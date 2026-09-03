"""Compare strict and guarded operator coverage on selected route roots."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from html import escape
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .chemistry import digest
from .generic_library import load_generic_library
from .generic_models import GenericTemplateLibrary
from .hierarchical_ranking import build_completion_prior_index
from .html_report import molecule_svg, reaction_svg
from .route_action_evaluation import (
    RouteActionEvaluationConfig,
    RouteActionSearcher,
    RouteActionStepEvaluation,
    evaluate_route_actions,
)
from .route_contract import ReactionRouteTree


OPERATOR_COVERAGE_COMPARISON_SCHEMA_VERSION = "1.0"
OPERATOR_COVERAGE_COMPARISON_ALGORITHM_VERSION = (
    "root_operator_coverage_comparison.v1"
)


@dataclass(frozen=True)
class OperatorLibrarySnapshot:
    """Identity and coverage denominators for one compared library."""

    label: str
    source_path: str
    source_sha256: str
    admission_policy: str
    source_row_count: int
    accepted_observation_count: int
    template_count: int
    operator_count: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible snapshot."""

        return asdict(self)


@dataclass(frozen=True)
class RootCoverageSummary:
    """Aggregate recovery counts for one library."""

    target_count: int
    candidate_coverage_count: int
    exact_recovery_count: int
    exact_top1_count: int
    strategy_recovery_count: int
    site_recovery_count: int
    operator_recovery_count: int
    any_source_patent_overlap_count: int
    exact_with_source_patent_overlap_count: int
    exact_without_source_patent_overlap_count: int

    def to_dict(self) -> dict[str, int]:
        """Return a JSON-compatible summary."""

        return asdict(self)


@dataclass(frozen=True)
class RootOperatorCoverageCase:
    """Strict and guarded replay results for one observed route root."""

    case_id: str
    selection_rank: int
    tree_id: str
    source_route_id: str | None
    patent_id: str | None
    split: str | None
    target_smiles: str
    observed_reaction_smiles: str
    observed_precursor_smiles: str | None
    strict: RouteActionStepEvaluation
    guarded: RouteActionStepEvaluation

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible comparison case."""

        return {
            "case_id": self.case_id,
            "selection_rank": self.selection_rank,
            "tree_id": self.tree_id,
            "source_route_id": self.source_route_id,
            "patent_id": self.patent_id,
            "split": self.split,
            "target_smiles": self.target_smiles,
            "observed_reaction_smiles": self.observed_reaction_smiles,
            "observed_precursor_smiles": self.observed_precursor_smiles,
            "strict": self.strict.to_dict(),
            "guarded": self.guarded.to_dict(),
        }


@dataclass(frozen=True)
class OperatorCoverageComparison:
    """Deterministic strict-versus-guarded root-action evaluation."""

    comparison_id: str
    source_review_path: str
    source_review_sha256: str
    strict_library: OperatorLibrarySnapshot
    guarded_library: OperatorLibrarySnapshot
    search_config: dict[str, Any]
    strict_summary: RootCoverageSummary
    guarded_summary: RootCoverageSummary
    cases: tuple[RootOperatorCoverageCase, ...]
    warnings: tuple[str, ...]
    schema_version: str = OPERATOR_COVERAGE_COMPARISON_SCHEMA_VERSION
    algorithm_version: str = OPERATOR_COVERAGE_COMPARISON_ALGORITHM_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible comparison."""

        return {
            "comparison_id": self.comparison_id,
            "source_review_path": self.source_review_path,
            "source_review_sha256": self.source_review_sha256,
            "strict_library": self.strict_library.to_dict(),
            "guarded_library": self.guarded_library.to_dict(),
            "search_config": self.search_config,
            "strict_summary": self.strict_summary.to_dict(),
            "guarded_summary": self.guarded_summary.to_dict(),
            "cases": [case.to_dict() for case in self.cases],
            "warnings": list(self.warnings),
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
        }


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _snapshot(
    label: str,
    path: Path,
    library: GenericTemplateLibrary,
) -> OperatorLibrarySnapshot:
    policy = str(library.definition.get("core_admission_policy") or "pass_only")
    return OperatorLibrarySnapshot(
        label=label,
        source_path=str(path.resolve()),
        source_sha256=_sha256(path),
        admission_policy=policy,
        source_row_count=library.source_row_count,
        accepted_observation_count=library.accepted_observation_count,
        template_count=len(library.templates),
        operator_count=len(library.operators),
    )


def _root_only_tree(tree: ReactionRouteTree) -> ReactionRouteTree:
    if tree.root.reaction is None:
        raise ValueError("operator coverage case requires a root reaction")
    children = tuple(
        replace(
            child,
            terminal=True,
            terminal_evidence="root_action_evaluation_frontier",
            unresolved_reason=None,
            reaction=None,
        )
        for child in tree.root.reaction.children
    )
    reaction = replace(tree.root.reaction, children=children)
    return replace(
        tree,
        root=replace(tree.root, reaction=reaction),
        reaction_count=1,
        maximum_depth=1,
    )


def _exact_overlap(step: RouteActionStepEvaluation) -> bool:
    return any(
        candidate.exact_precursor_match
        and candidate.source_patent_precedent_overlap
        for candidate in step.candidates
    )


def _summary(
    steps: tuple[RouteActionStepEvaluation, ...],
) -> RootCoverageSummary:
    exact_overlap = sum(_exact_overlap(step) for step in steps)
    exact = sum(step.exact_precursor_rank is not None for step in steps)
    return RootCoverageSummary(
        target_count=len(steps),
        candidate_coverage_count=sum(step.candidate_count > 0 for step in steps),
        exact_recovery_count=exact,
        exact_top1_count=sum(step.exact_precursor_rank == 1 for step in steps),
        strategy_recovery_count=sum(step.strategy_rank is not None for step in steps),
        site_recovery_count=sum(step.site_rank is not None for step in steps),
        operator_recovery_count=sum(step.operator_rank is not None for step in steps),
        any_source_patent_overlap_count=sum(
            step.source_patent_precedent_overlap for step in steps
        ),
        exact_with_source_patent_overlap_count=exact_overlap,
        exact_without_source_patent_overlap_count=exact - exact_overlap,
    )


def compare_root_operator_coverage(
    review: Mapping[str, Any],
    strict_library: GenericTemplateLibrary,
    guarded_library: GenericTemplateLibrary,
    *,
    source_review_path: str = "<in-memory>",
    source_review_sha256: str = "",
    strict_snapshot: OperatorLibrarySnapshot | None = None,
    guarded_snapshot: OperatorLibrarySnapshot | None = None,
    config: RouteActionEvaluationConfig = RouteActionEvaluationConfig(),
    searcher: RouteActionSearcher | None = None,
) -> OperatorCoverageComparison:
    """Replay selected route roots through strict and guarded libraries."""

    raw_cases = review.get("cases")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError("partition review must contain selected cases")
    strict_prior = (
        build_completion_prior_index(strict_library)
        if config.run_search and config.use_hierarchical_ranking
        else None
    )
    guarded_prior = (
        build_completion_prior_index(guarded_library)
        if config.run_search and config.use_hierarchical_ranking
        else None
    )
    cases = []
    for raw_case in raw_cases:
        raw_tree = raw_case.get("route_tree")
        if not isinstance(raw_tree, dict):
            raise ValueError("partition review case omits its route tree")
        tree = ReactionRouteTree.from_dict(raw_tree)
        root_tree = _root_only_tree(tree)
        evaluation_kwargs: dict[str, Any] = {
            "config": config,
        }
        if searcher is not None:
            evaluation_kwargs["searcher"] = searcher
        strict = evaluate_route_actions(
            root_tree,
            strict_library,
            completion_prior_index=strict_prior,
            **evaluation_kwargs,
        ).steps[0]
        guarded = evaluate_route_actions(
            root_tree,
            guarded_library,
            completion_prior_index=guarded_prior,
            **evaluation_kwargs,
        ).steps[0]
        cases.append(
            RootOperatorCoverageCase(
                case_id=str(raw_case.get("case_id") or tree.tree_id),
                selection_rank=int(raw_case.get("selection_rank") or len(cases) + 1),
                tree_id=tree.tree_id,
                source_route_id=tree.source_route_id,
                patent_id=tree.patent_id,
                split=tree.split,
                target_smiles=tree.target_smiles,
                observed_reaction_smiles=tree.root.reaction.reaction_smiles,
                observed_precursor_smiles=(
                    strict.observed_action.expected_precursor_smiles
                    if strict.observed_action.exact_precursors_verified
                    else None
                ),
                strict=strict,
                guarded=guarded,
            )
        )
    ordered = tuple(sorted(cases, key=lambda item: item.selection_rank))
    strict_steps = tuple(case.strict for case in ordered)
    guarded_steps = tuple(case.guarded for case in ordered)
    warnings = {
        "ROOT_ACTION_RECOVERY_IS_NOT_COMPLETE_ROUTE_VALIDATION",
        "VALIDATED_DEPARTURES_REMAINS_EXPERIMENTAL",
    }
    if any(case.split == "train" for case in ordered):
        warnings.add("SAMPLE_CONTAINS_TRAIN_ROUTES_NOT_BLIND_ACCURACY")
    if any(step.source_patent_precedent_overlap for step in guarded_steps):
        warnings.add("SOURCE_PATENT_PRECEDENT_OVERLAP_REPORTED_NOT_EXCLUDED")
    strict_identity = strict_snapshot or OperatorLibrarySnapshot(
        label="strict",
        source_path="<in-memory>",
        source_sha256="",
        admission_policy=str(
            strict_library.definition.get("core_admission_policy") or "pass_only"
        ),
        source_row_count=strict_library.source_row_count,
        accepted_observation_count=strict_library.accepted_observation_count,
        template_count=len(strict_library.templates),
        operator_count=len(strict_library.operators),
    )
    guarded_identity = guarded_snapshot or OperatorLibrarySnapshot(
        label="guarded",
        source_path="<in-memory>",
        source_sha256="",
        admission_policy=str(
            guarded_library.definition.get("core_admission_policy")
            or "pass_only"
        ),
        source_row_count=guarded_library.source_row_count,
        accepted_observation_count=guarded_library.accepted_observation_count,
        template_count=len(guarded_library.templates),
        operator_count=len(guarded_library.operators),
    )
    comparison_id = digest(
        "OPCOV1",
        source_review_sha256 or source_review_path,
        strict_identity.source_sha256 or strict_identity.admission_policy,
        guarded_identity.source_sha256 or guarded_identity.admission_policy,
        config.config_id,
        *(case.case_id for case in ordered),
    )
    return OperatorCoverageComparison(
        comparison_id=comparison_id,
        source_review_path=source_review_path,
        source_review_sha256=source_review_sha256,
        strict_library=strict_identity,
        guarded_library=guarded_identity,
        search_config=config.to_dict(),
        strict_summary=_summary(strict_steps),
        guarded_summary=_summary(guarded_steps),
        cases=ordered,
        warnings=tuple(sorted(warnings)),
    )


def build_root_operator_coverage_comparison(
    source_review: str | Path,
    strict_library_path: str | Path,
    guarded_library_path: str | Path,
    *,
    config: RouteActionEvaluationConfig = RouteActionEvaluationConfig(),
) -> OperatorCoverageComparison:
    """Load artifacts and compare strict and guarded root-action recovery."""

    review_path = Path(source_review)
    strict_path = Path(strict_library_path)
    guarded_path = Path(guarded_library_path)
    review = json.loads(review_path.read_text("utf-8"))
    if not isinstance(review, dict):
        raise ValueError("partition review must contain an object")
    strict_library = load_generic_library(strict_path)
    guarded_library = load_generic_library(guarded_path)
    return compare_root_operator_coverage(
        review,
        strict_library,
        guarded_library,
        source_review_path=str(review_path.resolve()),
        source_review_sha256=_sha256(review_path),
        strict_snapshot=_snapshot("strict pass-only", strict_path, strict_library),
        guarded_snapshot=_snapshot(
            "guarded validated departures",
            guarded_path,
            guarded_library,
        ),
        config=config,
    )


def _metric_table(comparison: OperatorCoverageComparison) -> str:
    rows = (
        ("Targets with candidates", "candidate_coverage_count"),
        ("Exact observed precursors", "exact_recovery_count"),
        ("Exact at rank 1", "exact_top1_count"),
        ("Observed strategy", "strategy_recovery_count"),
        ("Observed site", "site_recovery_count"),
        ("Observed operator", "operator_recovery_count"),
        ("Any source-patent overlap", "any_source_patent_overlap_count"),
        ("Exact without visible overlap", "exact_without_source_patent_overlap_count"),
    )
    return "".join(
        "<tr>"
        f"<th>{escape(label)}</th>"
        f"<td>{getattr(comparison.strict_summary, field)} / "
        f"{comparison.strict_summary.target_count}</td>"
        f"<td>{getattr(comparison.guarded_summary, field)} / "
        f"{comparison.guarded_summary.target_count}</td>"
        "</tr>"
        for label, field in rows
    )


def _rank(value: int | None) -> str:
    return "not recovered" if value is None else f"rank {value}"


def _candidate_cards(
    case: RootOperatorCoverageCase,
    step: RouteActionStepEvaluation,
) -> str:
    candidates = list(step.candidates[:5])
    exact = next(
        (item for item in step.candidates if item.exact_precursor_match),
        None,
    )
    if exact is not None and exact not in candidates:
        candidates.append(exact)
    if not candidates:
        return "<p class='empty'>No graph-validated candidates returned.</p>"
    values = []
    for candidate in candidates:
        tags = [f"rank {candidate.candidate_rank}", candidate.abstraction_level]
        if candidate.exact_precursor_match:
            tags.append("exact observed precursors")
        elif candidate.strategy_match:
            tags.append("observed strategy")
        elif candidate.site_match:
            tags.append("observed site")
        if candidate.source_patent_precedent_overlap:
            tags.append("source-patent overlap")
        values.append(
            "<article class='candidate'>"
            f"<div class='tags'>{''.join(f'<span>{escape(tag)}</span>' for tag in tags)}</div>"
            f"<div class='reaction'>{reaction_svg(f'{candidate.precursor_smiles}>>{case.target_smiles}')}</div>"
            f"<code>{escape(candidate.precursor_smiles)}</code>"
            f"<small>score {candidate.score:.3f} · "
            f"{candidate.independent_reference_support} independent references</small>"
            "</article>"
        )
    return "".join(values)


def _result_panel(
    case: RootOperatorCoverageCase,
    step: RouteActionStepEvaluation,
    label: str,
) -> str:
    return (
        "<section class='result-panel'>"
        f"<h3>{escape(label)}</h3>"
        f"<p><strong>{escape(step.outcome.replace('_', ' '))}</strong> · "
        f"{step.candidate_count} candidates</p>"
        "<div class='ranks'>"
        f"<span>exact: {_rank(step.exact_precursor_rank)}</span>"
        f"<span>strategy: {_rank(step.strategy_rank)}</span>"
        f"<span>site: {_rank(step.site_rank)}</span>"
        f"<span>operator: {_rank(step.operator_rank)}</span>"
        "</div>"
        f"<div class='candidate-list'>{_candidate_cards(case, step)}</div>"
        "</section>"
    )


def _case_panel(case: RootOperatorCoverageCase) -> str:
    return (
        f"<article class='case' data-split='{escape(case.split or 'unknown')}' "
        f"data-guarded='{escape(case.guarded.outcome)}'>"
        "<div class='case-heading'>"
        f"<div><span>Case {case.selection_rank}</span><h2>{escape(case.split or 'unknown')} route</h2>"
        f"<p>{escape(case.patent_id or 'unknown patent')} · "
        f"{escape(case.source_route_id or case.tree_id)}</p></div>"
        f"<strong>{escape(case.guarded.outcome.replace('_', ' '))}</strong></div>"
        "<div class='observed'>"
        f"<div>{molecule_svg(case.target_smiles, width=380, height=230)}</div>"
        f"<div><h3>Target and hidden observed answer</h3><code>{escape(case.target_smiles)}</code>"
        f"<details><summary>Reveal observed first reaction</summary>{reaction_svg(case.observed_reaction_smiles)}"
        f"<code>{escape(case.observed_reaction_smiles)}</code></details></div></div>"
        "<div class='comparison'>"
        f"{_result_panel(case, case.strict, 'Strict pass-only library')}"
        f"{_result_panel(case, case.guarded, 'Guarded validated-departures library')}"
        "</div><div class='review'><label>Guarded first step chemically reasonable?"
        f"<select data-case='{escape(case.case_id)}' data-field='decision'>"
        "<option value=''>Select…</option><option>yes</option><option>partly</option>"
        "<option>no</option><option>cannot_assess</option></select></label>"
        f"<label>Notes<textarea data-case='{escape(case.case_id)}' data-field='notes'></textarea></label>"
        "</div></article>"
    )


def render_operator_coverage_comparison_html(
    comparison: OperatorCoverageComparison,
) -> str:
    """Render a self-contained graphical strict-versus-guarded review."""

    payload = json.dumps(
        comparison.to_dict(), sort_keys=True, separators=(",", ":")
    ).replace("</", "<\\/")
    warning_items = "".join(
        f"<li>{escape(item)}</li>" for item in comparison.warnings
    )
    storage_key = f"operator-coverage:{comparison.comparison_id}"
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Root operator coverage comparison</title><style>
:root{{--ink:#203037;--muted:#697981;--paper:#fff;--bg:#edf2f3;--line:#d3dfe3;--blue:#28668a;--green:#267656;--amber:#9a650d}}*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,sans-serif}}header{{background:linear-gradient(125deg,#16384d,#286957);color:#fff;padding:28px max(20px,calc((100% - 1500px)/2))}}header p{{max-width:1050px;color:#dfecee}}main{{max-width:1500px;margin:auto;padding:18px}}code{{display:block;overflow-wrap:anywhere;font-size:11px}}table{{border-collapse:collapse;width:100%;background:#ffffff12}}th,td{{padding:7px 10px;border:1px solid #ffffff32;text-align:left}}.warnings{{background:#fff3da;color:#79500b;border:1px solid #efd79d;border-radius:9px;padding:8px 14px}}.toolbar{{position:sticky;top:0;z-index:4;display:flex;gap:12px;align-items:end;background:#fff;border:1px solid var(--line);border-radius:9px;padding:10px;margin:14px 0}}select,textarea,button{{font:inherit;border:1px solid #aebdc4;border-radius:6px;padding:7px;background:#fff}}button{{background:#245f79;color:#fff;border:0}}#visible{{margin-left:auto}}.case{{background:var(--paper);border:1px solid var(--line);border-radius:13px;padding:16px;margin:18px 0;box-shadow:0 4px 16px #18384410}}.case-heading{{display:flex;justify-content:space-between;gap:12px}}.case-heading span{{color:var(--blue);font-weight:800;text-transform:uppercase}}.case-heading h2{{margin:.2rem 0}}.observed{{display:grid;grid-template-columns:400px 1fr;gap:14px;align-items:center;border:1px solid var(--line);border-radius:9px;padding:10px;margin:12px 0}}.observed svg{{max-width:100%;height:auto;max-height:240px}}details summary{{cursor:pointer;font-weight:700}}.comparison{{display:grid;grid-template-columns:1fr 1fr;gap:12px}}.result-panel{{border:1px solid var(--line);border-top:5px solid var(--blue);border-radius:9px;padding:10px;min-width:0}}.result-panel h3{{margin:.1rem 0}}.ranks,.tags{{display:flex;flex-wrap:wrap;gap:5px}}.ranks span,.tags span{{background:#e7f0f4;border-radius:999px;padding:3px 7px;font-size:11px}}.candidate{{border-top:1px solid var(--line);padding:9px 0}}.candidate .reaction{{overflow:auto}}.candidate svg{{width:100%;height:auto;max-height:190px}}.candidate small{{display:block;color:var(--muted)}}.empty{{color:var(--muted)}}.review{{display:grid;grid-template-columns:280px 1fr;gap:10px;border-top:1px solid var(--line);margin-top:12px;padding-top:10px}}.review label{{display:grid;gap:4px;font-weight:700}}textarea{{min-height:70px}}@media(max-width:950px){{.comparison,.observed,.review{{grid-template-columns:1fr}}.toolbar{{position:static;flex-wrap:wrap}}#visible{{width:100%;margin:0}}}}
</style></head><body><header><p>EXPERIMENTAL OPERATOR ADMISSION · ROOT-ACTION REPLAY</p><h1>Strict versus guarded operator coverage</h1><p>Each target is searched independently. Returned candidates pass the same graph and forward-signature checks. The observed patent reaction is retained as a hidden reference answer; exact recovery is not the only chemically reasonable outcome and first-step recovery is not a complete route.</p><table><thead><tr><th>Measure</th><th>Strict</th><th>Guarded</th></tr></thead><tbody>{_metric_table(comparison)}</tbody></table></header><main><div class="warnings"><ul>{warning_items}</ul></div><div class="toolbar"><label>Split<select id="split"><option value="all">All</option><option value="train">Train</option><option value="validation">Validation</option><option value="test">Test</option></select></label><label>Guarded outcome<select id="outcome"><option value="all">All</option><option value="exact_top1">Exact top 1</option><option value="exact_lower_rank">Exact lower rank</option><option value="other_valid_candidates">Other candidates</option><option value="no_candidates">No candidates</option></select></label><button id="export">Export review JSON</button><span id="visible"></span></div>{''.join(_case_panel(case) for case in comparison.cases)}</main><script type="application/json" id="operator-coverage-comparison-data">{payload}</script><script>
const key={json.dumps(storage_key)},cards=[...document.querySelectorAll('.case')];let saved={{}};try{{saved=JSON.parse(localStorage.getItem(key)||'{{}}')}}catch(_){{saved={{}}}};for(const el of document.querySelectorAll('[data-field]')){{const id=el.dataset.case,field=el.dataset.field;el.value=(saved[id]||{{}})[field]||'';el.addEventListener(el.tagName==='TEXTAREA'?'input':'change',()=>{{saved[id]=saved[id]||{{}};saved[id][field]=el.value;try{{localStorage.setItem(key,JSON.stringify(saved))}}catch(_){{}}}})}}function filter(){{const split=document.querySelector('#split').value,outcome=document.querySelector('#outcome').value;let n=0;for(const card of cards){{const show=(split==='all'||card.dataset.split===split)&&(outcome==='all'||card.dataset.guarded===outcome);card.hidden=!show;if(show)n++}}document.querySelector('#visible').textContent=`${{n}} / ${{cards.length}} cases`}}document.querySelector('#split').addEventListener('change',filter);document.querySelector('#outcome').addEventListener('change',filter);document.querySelector('#export').addEventListener('click',()=>{{const out={{comparison_id:{json.dumps(comparison.comparison_id)},reviews:saved}};const blob=new Blob([JSON.stringify(out,null,2)],{{type:'application/json'}}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='operator-coverage-review.json';a.click();URL.revokeObjectURL(a.href)}});filter();
</script></body></html>"""


def write_operator_coverage_comparison(
    comparison: OperatorCoverageComparison,
    output_json: str | Path,
    output_html: str | Path,
) -> None:
    """Write deterministic JSON and graphical review artifacts."""

    json_path = Path(output_json)
    html_path = Path(output_html)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(comparison.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    html_path.write_text(
        render_operator_coverage_comparison_html(comparison),
        encoding="utf-8",
        newline="\n",
    )


__all__ = [
    "OPERATOR_COVERAGE_COMPARISON_ALGORITHM_VERSION",
    "OPERATOR_COVERAGE_COMPARISON_SCHEMA_VERSION",
    "OperatorCoverageComparison",
    "OperatorLibrarySnapshot",
    "RootCoverageSummary",
    "RootOperatorCoverageCase",
    "build_root_operator_coverage_comparison",
    "compare_root_operator_coverage",
    "render_operator_coverage_comparison_html",
    "write_operator_coverage_comparison",
]
