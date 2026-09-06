"""Baseline-preserving portfolio for validated latent-state routes.

The ordinary bounded planner is run without route-state guidance and retained
verbatim.  A second bounded search may contribute a small, explicitly separate
lane only when a novel route contains an operator supported by train-split
state-release evidence.  The exploratory lane cannot change baseline ranking.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import html
import json
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from .chemistry import digest
from .generic_models import GenericTemplateLibrary
from .multistep import (
    MultistepRetrosynthesisResult,
    MultistepRetrosynthesisRoute,
    plan_multistep_routes,
)
from .multistep_panel_review import (
    MultistepPanelCase,
    MultistepPanelTarget,
    render_multistep_panel_html,
)
from .html_report import molecule_svg, reaction_svg
from .route_state_learning import (
    LiteratureRouteActionSelector,
    LiteratureRouteOrderingGuidance,
    RouteStateLearningCatalog,
)


LATENT_ROUTE_PORTFOLIO_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "latent_route_portfolio.v1.json"
)
LATENT_ROUTE_PORTFOLIO_VERSION = "latent_route_portfolio.v1"


@dataclass(frozen=True)
class LatentRoutePortfolioPolicy:
    """Validated bounds for the separate exploratory route lane."""

    definition_id: str
    schema_version: str
    maximum_exploratory_routes: int
    require_supported_state_operator: bool
    require_route_novelty: bool
    preserve_baseline_routes: bool
    ranking_influence: str


@dataclass(frozen=True)
class LatentRouteAlternative:
    """One novel route carrying train-supported latent-state evidence."""

    route: MultistepRetrosynthesisRoute
    supported_state_operator_ids: tuple[str, ...]
    source_collection: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "route": self.route.to_dict(),
            "supported_state_operator_ids": list(
                self.supported_state_operator_ids
            ),
            "source_collection": self.source_collection,
        }


@dataclass(frozen=True)
class BaselinePreservingLatentPortfolio:
    """Ordinary result plus non-competing latent-state alternatives."""

    portfolio_id: str
    target_smiles: str
    baseline: MultistepRetrosynthesisResult
    exploratory_search: MultistepRetrosynthesisResult
    alternatives: tuple[LatentRouteAlternative, ...]
    route_state_definition_id: str
    portfolio_definition_id: str = LATENT_ROUTE_PORTFOLIO_VERSION
    schema_version: str = "1.0"

    @property
    def baseline_route_ids(self) -> tuple[str, ...]:
        return tuple(
            route.route_id
            for route in (self.baseline.routes or self.baseline.partial_routes)
        )

    @property
    def augmented_route_ids(self) -> tuple[str, ...]:
        return self.baseline_route_ids + tuple(
            item.route.route_id for item in self.alternatives
        )

    @property
    def baseline_preserved(self) -> bool:
        return self.augmented_route_ids[: len(self.baseline_route_ids)] == (
            self.baseline_route_ids
        )

    def alternative_result(self) -> MultistepRetrosynthesisResult:
        """Return an isolated result view for graphical comparison."""

        solved = tuple(item.route for item in self.alternatives if item.route.solved)
        partial = tuple(
            item.route for item in self.alternatives if not item.route.solved
        )
        return replace(
            self.exploratory_search,
            routes=solved,
            partial_routes=partial,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "portfolio_id": self.portfolio_id,
            "target_smiles": self.target_smiles,
            "baseline": self.baseline.to_dict(),
            "exploratory_search": self.exploratory_search.to_dict(),
            "alternatives": [item.to_dict() for item in self.alternatives],
            "baseline_route_ids": list(self.baseline_route_ids),
            "augmented_route_ids": list(self.augmented_route_ids),
            "baseline_preserved": self.baseline_preserved,
            "route_state_definition_id": self.route_state_definition_id,
            "portfolio_definition_id": self.portfolio_definition_id,
            "schema_version": self.schema_version,
            "warnings": [
                "EXPLORATORY_ROUTES_DO_NOT_CHANGE_BASELINE_RANKING",
                "STATE_OPERATOR_SUPPORT_IS_NOT_PROOF_THAT_PROTECTION_IS_NEEDED",
                "EVERY_RETAINED_ACTION_PASSED_THE_CANONICAL_PLANNER_GATES",
            ],
        }


def load_latent_route_portfolio_policy(
    source: str | Path = LATENT_ROUTE_PORTFOLIO_DEFINITION_PATH,
) -> LatentRoutePortfolioPolicy:
    """Load and validate the baseline-preserving portfolio policy."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if value.get("definition_id") != LATENT_ROUTE_PORTFOLIO_VERSION:
        raise ValueError("unexpected latent-route portfolio definition ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported latent-route portfolio schema")
    maximum = value.get("maximum_exploratory_routes")
    if isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1:
        raise ValueError("maximum exploratory routes must be positive")
    for key in (
        "require_supported_state_operator",
        "require_route_novelty",
        "preserve_baseline_routes",
    ):
        if value.get(key) is not True:
            raise ValueError(f"{key} must remain enabled")
    if value.get("ranking_influence") != "separate_opt_in_exploratory_lane":
        raise ValueError("latent alternatives cannot modify baseline ranking")
    return LatentRoutePortfolioPolicy(**value)


def _route_identity(route: MultistepRetrosynthesisRoute) -> tuple[str, ...]:
    return tuple(step.candidate.proposed_reaction_smiles for step in route.steps)


def _supported_route_operators(
    route: MultistepRetrosynthesisRoute,
    supported: frozenset[str],
) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                step.candidate.operator_id
                for step in route.steps
                if step.candidate.operator_id in supported
            }
        )
    )


def build_baseline_preserving_latent_portfolio(
    baseline: MultistepRetrosynthesisResult,
    exploratory_search: MultistepRetrosynthesisResult,
    catalog: RouteStateLearningCatalog,
    *,
    policy: LatentRoutePortfolioPolicy | None = None,
) -> BaselinePreservingLatentPortfolio:
    """Select novel supported alternatives without changing the baseline."""

    active_policy = policy or load_latent_route_portfolio_policy()
    if baseline.target_smiles != exploratory_search.target_smiles:
        raise ValueError("baseline and exploratory targets must match")
    baseline_identities = {
        _route_identity(route)
        for route in baseline.routes + baseline.partial_routes
    }
    alternatives: list[LatentRouteAlternative] = []
    supported = catalog.supported_state_operator_ids
    for source_collection, routes in (
        ("solved", exploratory_search.routes),
        ("partial", exploratory_search.partial_routes),
    ):
        for route in routes:
            if _route_identity(route) in baseline_identities:
                continue
            operators = _supported_route_operators(route, supported)
            if not operators:
                continue
            alternatives.append(
                LatentRouteAlternative(
                    route=route,
                    supported_state_operator_ids=operators,
                    source_collection=source_collection,
                )
            )
            if len(alternatives) >= active_policy.maximum_exploratory_routes:
                break
        if len(alternatives) >= active_policy.maximum_exploratory_routes:
            break
    portfolio_id = digest(
        "LRP1",
        baseline.target_smiles,
        catalog.definition_id,
        *(
            route.route_id
            for route in baseline.routes + baseline.partial_routes
        ),
        *(item.route.route_id for item in alternatives),
    )
    result = BaselinePreservingLatentPortfolio(
        portfolio_id=portfolio_id,
        target_smiles=baseline.target_smiles,
        baseline=baseline,
        exploratory_search=exploratory_search,
        alternatives=tuple(alternatives),
        route_state_definition_id=catalog.definition_id,
        portfolio_definition_id=active_policy.definition_id,
    )
    if not result.baseline_preserved:
        raise AssertionError("baseline route order was not preserved")
    return result


Planner = Callable[..., MultistepRetrosynthesisResult]


def plan_baseline_preserving_latent_portfolio(
    target_smiles: str,
    library: GenericTemplateLibrary,
    literature_index: Any,
    catalog: RouteStateLearningCatalog,
    *,
    enable_experimental_ordering: bool = False,
    planner: Planner = plan_multistep_routes,
    **planner_kwargs: Any,
) -> BaselinePreservingLatentPortfolio:
    """Run ordinary and latent searches independently, then keep novel lanes."""

    forbidden = {"route_action_selector", "search_guidance"}.intersection(
        planner_kwargs
    )
    if forbidden:
        raise ValueError(
            "portfolio owns planner hooks: " + ", ".join(sorted(forbidden))
        )
    baseline = planner(
        target_smiles,
        library,
        literature_index,
        **planner_kwargs,
    )
    exploratory = planner(
        target_smiles,
        library,
        literature_index,
        route_action_selector=LiteratureRouteActionSelector(catalog),
        search_guidance=(
            LiteratureRouteOrderingGuidance(catalog)
            if enable_experimental_ordering
            else None
        ),
        **planner_kwargs,
    )
    return build_baseline_preserving_latent_portfolio(
        baseline, exploratory, catalog
    )


def write_latent_route_portfolio(
    portfolio: BaselinePreservingLatentPortfolio,
    output_json: str | Path,
    output_html: str | Path,
    *,
    title: str = "Baseline-preserving latent-route portfolio",
) -> dict[str, Any]:
    """Write the typed portfolio and a structure-rich comparison."""

    json_path = Path(output_json)
    html_path = Path(output_html)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(portfolio.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    panel = MultistepPanelCase(
        target=MultistepPanelTarget(
            target_id=portfolio.portfolio_id,
            name="Target",
            category="latent_state_portfolio",
            smiles=portfolio.target_smiles,
        ),
        baseline=portfolio.baseline,
        policy=portfolio.alternative_result(),
        evaluation_metrics=(
            ("Baseline preserved", str(portfolio.baseline_preserved)),
            ("Novel latent routes", str(len(portfolio.alternatives))),
        ),
    )
    document = render_multistep_panel_html(
        (panel,),
        title=title,
        top_k=max(1, len(portfolio.baseline_route_ids)),
        metadata={
            "panel_id": portfolio.portfolio_id,
            "warnings": (
                "The latent-state lane is additional and never replaces a baseline route.",
                "Operator support does not establish that protection is necessary.",
            ),
        },
        policy_label="Novel latent-state lane",
    )
    html_path.write_text(document, encoding="utf-8", newline="\n")
    return {
        "portfolio_id": portfolio.portfolio_id,
        "baseline_preserved": portfolio.baseline_preserved,
        "baseline_route_count": len(portfolio.baseline_route_ids),
        "novel_latent_route_count": len(portfolio.alternatives),
        "output_json": str(json_path.resolve()),
        "output_html": str(html_path.resolve()),
    }


def _load_json_mapping(source: str | Path | Mapping[str, Any]) -> Mapping[str, Any]:
    if isinstance(source, Mapping):
        return source
    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError("evaluation artifact must contain a JSON object")
    return value


def _raw_selected_routes(result: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    routes = result.get("routes") or result.get("partial_routes") or ()
    return [item for item in routes if isinstance(item, Mapping)]


def _raw_route_identity(route: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(
        str((step.get("candidate") or {}).get("proposed_reaction_smiles") or "")
        for step in route.get("steps") or ()
        if isinstance(step, Mapping)
    )


def _raw_supported_operators(
    route: Mapping[str, Any], supported: frozenset[str]
) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                str((step.get("candidate") or {}).get("operator_id") or "")
                for step in route.get("steps") or ()
                if isinstance(step, Mapping)
                and str((step.get("candidate") or {}).get("operator_id") or "")
                in supported
            }
        )
    )


def build_dataset_latent_portfolio_review(
    baseline_evaluation: str | Path | Mapping[str, Any],
    exploratory_evaluation: str | Path | Mapping[str, Any],
    catalog: RouteStateLearningCatalog,
) -> dict[str, Any]:
    """Combine fixed-compute evaluation outputs without replacing baselines."""

    baseline = _load_json_mapping(baseline_evaluation)
    exploratory = _load_json_mapping(exploratory_evaluation)
    baseline_cases = baseline.get("cases") or ()
    exploratory_cases = exploratory.get("cases") or ()
    if len(baseline_cases) != len(exploratory_cases) or not baseline_cases:
        raise ValueError("evaluation artifacts must contain aligned cases")
    supported = catalog.supported_state_operator_ids
    cases = []
    novel_count = 0
    for baseline_case, exploratory_case in zip(
        baseline_cases, exploratory_cases, strict=True
    ):
        if baseline_case.get("case_id") != exploratory_case.get("case_id"):
            raise ValueError("evaluation case identities are not aligned")
        baseline_panel = baseline_case.get("panel_case") or {}
        exploratory_panel = exploratory_case.get("panel_case") or {}
        baseline_result = baseline_panel.get("baseline") or {}
        exploratory_result = exploratory_panel.get("baseline") or {}
        baseline_routes = _raw_selected_routes(baseline_result)
        identities = {_raw_route_identity(route) for route in baseline_routes}
        alternative = None
        alternative_operators: tuple[str, ...] = ()
        for route in _raw_selected_routes(exploratory_result):
            operators = _raw_supported_operators(route, supported)
            if operators and _raw_route_identity(route) not in identities:
                alternative = route
                alternative_operators = operators
                novel_count += 1
                break
        cases.append(
            {
                "case_id": baseline_case.get("case_id"),
                "selection_rank": baseline_case.get("selection_rank"),
                "target": baseline_panel.get("target") or {},
                "reference_reactions": list(
                    baseline_panel.get("reference_reactions") or ()
                ),
                "baseline_routes": baseline_routes,
                "baseline_route_ids": [
                    str(route.get("route_id") or "") for route in baseline_routes
                ],
                "latent_alternative": alternative,
                "supported_state_operator_ids": list(alternative_operators),
                "baseline_preserved": True,
            }
        )
    return {
        "schema_version": "1.0",
        "definition_id": LATENT_ROUTE_PORTFOLIO_VERSION,
        "review_id": digest(
            "LRPD1",
            str(baseline.get("evaluation_id") or ""),
            str(exploratory.get("evaluation_id") or ""),
            catalog.definition_id,
        ),
        "baseline_evaluation_id": baseline.get("evaluation_id"),
        "exploratory_evaluation_id": exploratory.get("evaluation_id"),
        "route_state_definition_id": catalog.definition_id,
        "summary": {
            "target_count": len(cases),
            "baseline_preserved_count": len(cases),
            "targets_with_novel_latent_route": novel_count,
            "baseline_metrics_lower_bound": dict(baseline.get("summary") or {}),
        },
        "cases": cases,
        "warnings": [
            "LATENT_ROUTES_ARE_ADDITIONAL_AND_DO_NOT_REPLACE_BASELINE_ROUTES",
            "USEFULNESS_REQUIRES_BLIND_CHEMIST_REVIEW",
            "STATE_OPERATOR_SUPPORT_DOES_NOT_PROVE_PROTECTION_IS_NEEDED",
        ],
    }


def _raw_route_card(
    route: Mapping[str, Any],
    *,
    label: str,
    supported: Sequence[str] = (),
) -> str:
    steps = []
    for index, step in enumerate(route.get("steps") or (), start=1):
        candidate = step.get("candidate") or {}
        reaction = str(candidate.get("proposed_reaction_smiles") or "")
        operator = str(candidate.get("operator_id") or "")
        marker = " · supported state operator" if operator in supported else ""
        steps.append(
            "<div class='reaction'>"
            f"{reaction_svg(reaction)}"
            f"<p>Step {index} · {html.escape(str(candidate.get('strategic_class') or 'unresolved'))}{marker}</p>"
            "</div>"
        )
    status = "solved" if route.get("solved") else "partial"
    return (
        "<article class='route'>"
        f"<h4>{html.escape(label)} <span>{status}</span></h4>"
        f"<p>{len(steps)} reaction(s) · cost {float(route.get('route_cost') or 0.0):.3f}</p>"
        f"{''.join(steps)}</article>"
    )


def render_dataset_latent_portfolio_html(
    review: Mapping[str, Any],
    *,
    title: str = "Dataset-10 baseline-preserving latent-route review",
) -> str:
    """Render baseline routes beside the additional latent-state lane."""

    cards = []
    for case in review.get("cases") or ():
        target = case.get("target") or {}
        target_smiles = str(target.get("smiles") or "")
        references = "".join(
            f"<div class='reaction'>{reaction_svg(str(reaction))}</div>"
            for reaction in case.get("reference_reactions") or ()
        )
        baseline_routes = "".join(
            _raw_route_card(route, label=f"Baseline route {index}")
            for index, route in enumerate(case.get("baseline_routes") or (), 1)
        )
        alternative = case.get("latent_alternative")
        latent = (
            _raw_route_card(
                alternative,
                label="Additional latent-state route",
                supported=case.get("supported_state_operator_ids") or (),
            )
            if isinstance(alternative, Mapping)
            else "<p class='empty'>No novel supported latent route in this bounded run.</p>"
        )
        cards.append(
            f"<section class='case' data-case-id='{html.escape(str(case.get('case_id') or ''))}'>"
            f"<h2>Case {html.escape(str(case.get('selection_rank') or ''))}</h2>"
            f"<div class='target'>{molecule_svg(target_smiles, width=400, height=210)}<code>{html.escape(target_smiles)}</code></div>"
            "<div class='columns'>"
            f"<div><h3>Observed precedent</h3>{references}</div>"
            f"<div><h3>Planner baseline — preserved</h3>{baseline_routes}</div>"
            f"<div class='latent'><h3>Novel latent-state lane</h3>{latent}</div>"
            "</div><div class='review'><label>Latent-route utility<select class='verdict'><option value='unreviewed'>Unreviewed</option><option value='useful'>Useful addition</option><option value='unnecessary'>Unnecessary protection</option><option value='questionable'>Questionable chemistry</option></select></label><label>Notes<textarea class='notes' rows='2'></textarea></label></div></section>"
        )
    summary = review.get("summary") or {}
    notices = "".join(
        f"<li>{html.escape(str(item))}</li>" for item in review.get("warnings") or ()
    )
    return f"""<!doctype html><html lang='en'><head><meta charset='utf-8'>
<meta name='viewport' content='width=device-width,initial-scale=1'><title>{html.escape(title)}</title>
<style>
body{{margin:0;background:#f4f6f5;color:#17221e;font:14px/1.45 system-ui,Segoe UI,sans-serif}}header{{background:#173f33;color:white;padding:24px}}main{{max-width:1580px;margin:auto;padding:18px}}.summary{{display:flex;gap:12px;flex-wrap:wrap}}.summary b{{background:#285646;padding:8px 12px;border-radius:7px}}button{{border:0;border-radius:6px;padding:8px 12px;background:#e8f2ed;color:#173f33;font-weight:700;cursor:pointer}}.case{{background:white;border:1px solid #d5dfda;border-radius:12px;margin:18px 0;padding:16px}}.target{{display:flex;align-items:center;gap:14px;border-bottom:1px solid #dbe3df}}.target svg{{height:180px;max-width:400px}}code{{overflow-wrap:anywhere}}.columns{{display:grid;grid-template-columns:repeat(3,minmax(360px,1fr));gap:14px;margin-top:12px}}.columns>div{{border:1px solid #d7dfdb;border-radius:9px;padding:10px;min-width:0}}.latent{{background:#f0f8f3;border-color:#9bc7ad!important}}.route{{border-top:1px solid #dce3df;padding-top:8px;margin-top:8px}}.route h4{{margin:0}}.route h4 span{{font-size:12px;background:#e8f0ec;border-radius:999px;padding:3px 7px}}.reaction{{overflow:auto;background:#fafbfa;border-radius:7px;margin:7px 0}}.reaction svg{{display:block;width:100%;height:145px}}.reaction p{{margin:0;padding:4px 8px;color:#53625b}}.empty{{color:#6b756f}}.review{{display:grid;grid-template-columns:220px 1fr;gap:12px;border-top:1px solid #dbe3df;padding-top:12px;margin-top:12px}}.review label{{display:grid;gap:4px;color:#53625b}}.review select,.review textarea{{border:1px solid #bcc9c2;border-radius:6px;padding:7px;background:white}}@media(max-width:1100px){{.columns{{grid-template-columns:1fr}}}}@media(max-width:700px){{.review{{grid-template-columns:1fr}}}}
</style></head><body><header><h1>{html.escape(title)}</h1><div class='summary'><b>{summary.get('target_count',0)} targets</b><b>{summary.get('baseline_preserved_count',0)} baselines preserved</b><b>{summary.get('targets_with_novel_latent_route',0)} novel latent lanes</b><button id='export' type='button'>Export review JSON</button></div><ul>{notices}</ul></header><main>{''.join(cards)}</main><script>
const key={json.dumps('latent-route-review-' + str(review.get('review_id') or 'v1'))};
const saved=JSON.parse(localStorage.getItem(key)||'{{}}');
const cards=[...document.querySelectorAll('.case')];
for(const card of cards){{const id=card.dataset.caseId;const value=saved[id]||{{}};const verdict=card.querySelector('.verdict');const notes=card.querySelector('.notes');verdict.value=value.verdict||'unreviewed';notes.value=value.notes||'';const persist=()=>{{saved[id]={{verdict:verdict.value,notes:notes.value}};localStorage.setItem(key,JSON.stringify(saved));}};verdict.addEventListener('change',persist);notes.addEventListener('input',persist);}}
document.querySelector('#export').addEventListener('click',()=>{{const payload=cards.map(card=>({{case_id:card.dataset.caseId,...(saved[card.dataset.caseId]||{{verdict:'unreviewed',notes:''}})}}));const blob=new Blob([JSON.stringify(payload,null,2)],{{type:'application/json'}});const link=document.createElement('a');link.href=URL.createObjectURL(blob);link.download='latent-route-review.json';link.click();URL.revokeObjectURL(link.href);}});
</script></body></html>"""


def write_dataset_latent_portfolio_review(
    review: Mapping[str, Any],
    output_json: str | Path,
    output_html: str | Path,
) -> dict[str, Any]:
    """Write a dataset-level baseline-plus-latent portfolio review."""

    json_path = Path(output_json)
    html_path = Path(output_html)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(review, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    html_path.write_text(
        render_dataset_latent_portfolio_html(review),
        encoding="utf-8",
        newline="\n",
    )
    return {
        "review_id": review.get("review_id"),
        "summary": review.get("summary"),
        "output_json": str(json_path.resolve()),
        "output_html": str(html_path.resolve()),
    }


__all__ = [
    "BaselinePreservingLatentPortfolio",
    "LATENT_ROUTE_PORTFOLIO_DEFINITION_PATH",
    "LATENT_ROUTE_PORTFOLIO_VERSION",
    "LatentRouteAlternative",
    "LatentRoutePortfolioPolicy",
    "build_baseline_preserving_latent_portfolio",
    "build_dataset_latent_portfolio_review",
    "load_latent_route_portfolio_policy",
    "plan_baseline_preserving_latent_portfolio",
    "render_dataset_latent_portfolio_html",
    "write_dataset_latent_portfolio_review",
    "write_latent_route_portfolio",
]
