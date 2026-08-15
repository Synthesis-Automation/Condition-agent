"""Chemist-facing comparison panel for bounded multistep route searches."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import html
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

from rdkit import Chem

from .html_report import molecule_svg, reaction_svg
from .multistep import MultistepRetrosynthesisResult, MultistepRetrosynthesisRoute


_DEFAULT_PANEL_PATH = (
    Path(__file__).resolve().parent
    / "definitions"
    / "multistep_familiar_target_panel.v1.json"
)


@dataclass(frozen=True)
class MultistepPanelTarget:
    """One named target in a qualitative multistep review panel."""

    target_id: str
    name: str
    category: str
    smiles: str

    def __post_init__(self) -> None:
        if not self.target_id or not self.name or not self.category:
            raise ValueError("multistep panel target metadata cannot be empty")
        molecule = Chem.MolFromSmiles(self.smiles)
        if molecule is None:
            raise ValueError(f"invalid multistep panel SMILES: {self.smiles}")


@dataclass(frozen=True)
class MultistepPanelCase:
    """Baseline and optional policy search results for one target."""

    target: MultistepPanelTarget
    baseline: MultistepRetrosynthesisResult
    policy: Optional[MultistepRetrosynthesisResult] = None
    reference_route_id: Optional[str] = None
    reference_reactions: tuple[str, ...] = ()
    reference_maximum_depth: Optional[int] = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible panel case."""

        return {
            "target": asdict(self.target),
            "baseline": self.baseline.to_dict(),
            "policy": self.policy.to_dict() if self.policy is not None else None,
            "reference_route_id": self.reference_route_id,
            "reference_reactions": list(self.reference_reactions),
            "reference_maximum_depth": self.reference_maximum_depth,
            "same_ranked_routes": (
                _route_projection(self.baseline) == _route_projection(self.policy)
                if self.policy is not None
                else None
            ),
        }


def load_multistep_panel_targets(
    source: str | Path = _DEFAULT_PANEL_PATH,
) -> tuple[MultistepPanelTarget, ...]:
    """Load and validate a versioned familiar-target panel."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported multistep target-panel schema")
    panel_id = str(value.get("panel_id") or "")
    if not panel_id.startswith(
        (
            "multistep_familiar_target_panel.v1@",
            "multistep_complex_route_panel.v1@",
            "multistep_route_policy_test_panel.v1@",
            "multistep_route_policy_validation_panel.v1@",
        )
    ):
        raise ValueError("unexpected multistep target-panel identity")
    targets = tuple(
        MultistepPanelTarget(
            target_id=str(item["target_id"]),
            name=str(item["name"]),
            category=str(item["category"]),
            smiles=str(item["smiles"]),
        )
        for item in value.get("targets") or ()
    )
    if not targets:
        raise ValueError("multistep target panel cannot be empty")
    identities = tuple(item.target_id for item in targets)
    if len(set(identities)) != len(identities):
        raise ValueError("multistep target IDs must be unique")
    return targets


def _selected_routes(
    result: Optional[MultistepRetrosynthesisResult],
) -> tuple[MultistepRetrosynthesisRoute, ...]:
    if result is None:
        return ()
    return result.routes or result.partial_routes


def _route_projection(
    result: Optional[MultistepRetrosynthesisResult],
) -> tuple[tuple[str, float], ...]:
    return tuple(
        (route.route_id, route.route_cost) for route in _selected_routes(result)
    )


def _route_card(route: MultistepRetrosynthesisRoute, rank: int) -> str:
    status = "Solved" if route.solved else "Partial"
    steps = "".join(
        '<div class="reaction">'
        + reaction_svg(step.candidate.proposed_reaction_smiles)
        + '<div class="step-meta">'
        + f"depth {step.depth + 1} · cost {step.step_cost:.3f}"
        + (
            f" · policy p={step.route_policy_probability:.3f}"
            if step.route_policy_probability is not None
            else ""
        )
        + "</div></div>"
        for step in route.steps
    )
    leaves = "".join(
        '<span class="leaf">'
        f"{html.escape(leaf.smiles)} · {html.escape(leaf.terminal_evidence)}"
        "</span>"
        for leaf in route.leaves
    )
    warnings = "".join(
        f"<li>{html.escape(warning)}</li>" for warning in route.warnings
    )
    return (
        '<article class="route">'
        f'<div class="route-head"><strong>Route {rank}</strong>'
        f'<span class="badge">{status}</span>'
        f'<span>{route.reaction_count} reaction(s)</span>'
        f'<span>cost {route.route_cost:.3f}</span></div>'
        f'<div class="route-steps">{steps}</div>'
        f'<div class="leaves"><b>Starting-material assessment:</b>{leaves}</div>'
        + (f'<ul class="warnings">{warnings}</ul>' if warnings else "")
        + "</article>"
    )


def _method_column(
    label: str,
    result: Optional[MultistepRetrosynthesisResult],
    *,
    top_k: int,
) -> str:
    if result is None:
        return f'<section class="method"><h3>{html.escape(label)}</h3><p>Not run.</p></section>'
    routes = _selected_routes(result)[:top_k]
    status = (
        f"{len(result.routes)} solved"
        if result.routes
        else f"0 solved · {len(result.partial_routes)} partial"
    )
    policy_note = ""
    if result.route_action_policy_model_id:
        policy_note = (
            '<p class="policy-state '
            + ("active" if result.route_action_policy_active else "inactive")
            + '">Policy '
            + ("active" if result.route_action_policy_active else "inactive")
            + f" · scored {result.diagnostics.route_policy_scored_actions} actions"
            + f" · reordered {result.diagnostics.route_policy_reordered_expansions}</p>"
        )
    cards = "".join(_route_card(route, rank) for rank, route in enumerate(routes, 1))
    if not cards:
        cards = '<p class="empty">No route or bounded partial route returned.</p>'
    return (
        f'<section class="method"><h3>{html.escape(label)}</h3>'
        f'<p class="method-summary">{status}</p>{policy_note}{cards}</section>'
    )


def _reference_column(case: MultistepPanelCase) -> str:
    if not case.reference_reactions:
        return ""
    reactions = "".join(
        '<div class="reaction">'
        + reaction_svg(reaction)
        + f'<div class="step-meta">recorded reaction {index}</div></div>'
        for index, reaction in enumerate(case.reference_reactions, 1)
    )
    depth = (
        f" · maximum depth {case.reference_maximum_depth}"
        if case.reference_maximum_depth is not None
        else ""
    )
    return (
        '<section class="method reference"><h3>Observed precedent</h3>'
        f'<p class="method-summary">{len(case.reference_reactions)} recorded reactions'
        f'{depth}</p><p class="reference-id">{html.escape(case.reference_route_id or "")}</p>'
        f'<div class="route-steps">{reactions}</div></section>'
    )


def _case_card(case: MultistepPanelCase, index: int, *, top_k: int) -> str:
    if case.policy is None:
        comparison = "Observed vs planner"
        comparison_class = "reference"
    else:
        same = _route_projection(case.baseline) == _route_projection(case.policy)
        comparison = "Same ranked routes" if same else "Ranking changed"
        comparison_class = "same" if same else "changed"
    search_text = " ".join(
        (
            case.target.target_id,
            case.target.name,
            case.target.category,
            case.target.smiles,
        )
    ).lower()
    return (
        f'<article class="case" data-target-id="{html.escape(case.target.target_id)}" '
        f'data-category="{html.escape(case.target.category)}" '
        f'data-search="{html.escape(search_text)}">'
        '<div class="target-head">'
        f'<div><span class="case-number">{index:02d}</span>'
        f'<h2>{html.escape(case.target.name)}</h2>'
        f'<span class="category">{html.escape(case.target.category)}</span></div>'
        f'<span class="comparison {comparison_class}">{comparison}</span></div>'
        '<div class="target-structure">'
        + molecule_svg(case.target.smiles, width=420, height=220)
        + f'<code>{html.escape(case.target.smiles)}</code></div>'
        '<div class="methods">'
        + _reference_column(case)
        + _method_column("Baseline", case.baseline, top_k=top_k)
        + (
            _method_column("Route policy", case.policy, top_k=top_k)
            if case.policy is not None
            else ""
        )
        + "</div>"
        '<div class="review"><label>Status<select class="review-status">'
        '<option value="unreviewed">Unreviewed</option><option value="accept">Accept</option>'
        '<option value="question">Questionable</option><option value="reject">Reject</option>'
        '</select></label><label>Notes<textarea class="review-note" rows="2" '
        'placeholder="Chemistry, ranking, missing disconnection…"></textarea></label></div>'
        "</article>"
    )


def _stylesheet() -> str:
    return """
:root{--ink:#18201d;--muted:#66726d;--line:#d8dfdb;--paper:#f4f6f4;--green:#176b4a;--amber:#925f00;--red:#a33b2f}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}
header{background:#173f33;color:white;padding:26px max(24px,calc((100% - 1500px)/2));display:flex;justify-content:space-between;gap:24px;align-items:end}
header h1{margin:0 0 5px;font-size:27px}header p{margin:0;color:#dcebe5;max-width:850px}.summary{display:grid;grid-template-columns:repeat(5,minmax(120px,1fr));gap:1px;background:#42675b}.summary div{background:#234e41;padding:11px 15px}.summary span{display:block;color:#bcd2ca;font-size:12px}.summary strong{font-size:20px}
main{max-width:1500px;margin:auto;padding:18px}.toolbar{position:sticky;top:0;z-index:5;background:#fff;border:1px solid var(--line);border-radius:10px;padding:10px 12px;display:flex;gap:12px;align-items:end;box-shadow:0 3px 14px #203c3120}.toolbar label{display:grid;gap:3px;color:var(--muted);font-size:12px}.toolbar input,.toolbar select,.review select,.review textarea{border:1px solid #bec9c4;border-radius:6px;padding:7px;background:#fff}.toolbar input{width:260px}.toolbar button{border:0;border-radius:6px;background:#315e50;color:white;padding:8px 12px;cursor:pointer}.visible{margin-left:auto;color:var(--muted)}
.case{background:white;border:1px solid var(--line);border-radius:12px;margin:18px 0;padding:16px;box-shadow:0 3px 13px #223b3210}.target-head{display:flex;justify-content:space-between;align-items:center;gap:12px}.target-head>div{display:flex;align-items:center;gap:10px}.target-head h2{margin:0;font-size:22px}.case-number{color:var(--muted);font-variant-numeric:tabular-nums}.category,.comparison,.badge{border-radius:999px;padding:3px 8px;font-size:12px;background:#edf2ef}.comparison.same{color:var(--green)}.comparison.changed{background:#fff0d2;color:var(--amber)}.comparison.reference{background:#eee9dc;color:#66583b}
.target-structure{display:flex;align-items:center;gap:18px;overflow:auto;border-bottom:1px solid var(--line);padding:8px 0}.target-structure svg{max-width:420px;height:190px}.target-structure code{color:#43514c;overflow-wrap:anywhere}.methods{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:16px;margin-top:14px}.method{min-width:0;border:1px solid var(--line);border-radius:9px;padding:12px}.method.reference{background:#f7f5ef}.method h3{margin:0;font-size:18px}.method-summary,.reference-id{color:var(--muted);margin:2px 0 9px}.reference-id{font-family:ui-monospace,monospace}.policy-state{border-left:4px solid var(--green);padding:5px 8px;background:#eff8f3}.policy-state.inactive{border-color:var(--amber);background:#fff8e8}.route{border-top:1px solid var(--line);padding-top:10px;margin-top:10px}.route-head{display:flex;gap:9px;align-items:center;color:var(--muted)}.route-head strong{color:var(--ink)}.route-steps{display:grid;gap:7px}.reaction{overflow:auto;background:#fafbfa;border-radius:7px}.reaction svg{width:100%;min-width:520px;height:auto;max-height:220px}.step-meta{padding:0 8px 6px;color:var(--muted);font-size:12px}.leaves{display:flex;gap:5px;flex-wrap:wrap;align-items:center;margin-top:7px;font-size:12px}.leaf{background:#edf4f0;border-radius:5px;padding:3px 6px;overflow-wrap:anywhere}.warnings{margin:7px 0;color:var(--amber)}
.review{display:grid;grid-template-columns:180px 1fr;gap:12px;border-top:1px solid var(--line);margin-top:14px;padding-top:12px}.review label{display:grid;gap:4px;color:var(--muted);font-size:12px}.review textarea{width:100%;resize:vertical}.case[data-review="accept"]{border-left:5px solid var(--green)}.case[data-review="question"]{border-left:5px solid var(--amber)}.case[data-review="reject"]{border-left:5px solid var(--red)}
@media(max-width:900px){.summary{grid-template-columns:repeat(2,1fr)}.methods{grid-template-columns:1fr}.toolbar{position:static;flex-wrap:wrap}.toolbar input{width:190px}.review{grid-template-columns:1fr}}
"""


def _script(storage_key: str) -> str:
    return f"""
const key={json.dumps(storage_key)};
const saved=JSON.parse(localStorage.getItem(key)||'{{}}');
const cards=[...document.querySelectorAll('.case')];
for(const card of cards){{
  const id=card.dataset.targetId, state=saved[id]||{{}};
  const status=card.querySelector('.review-status'), note=card.querySelector('.review-note');
  status.value=state.status||'unreviewed'; note.value=state.note||''; card.dataset.review=status.value;
  const persist=()=>{{saved[id]={{status:status.value,note:note.value}};card.dataset.review=status.value;localStorage.setItem(key,JSON.stringify(saved));filter();}};
  status.addEventListener('change',persist);note.addEventListener('input',persist);
}}
function filter(){{
 const query=document.querySelector('#search').value.toLowerCase();
 const category=document.querySelector('#category').value;
 const status=document.querySelector('#status').value;let visible=0;
 for(const card of cards){{const show=(!query||card.dataset.search.includes(query))&&(category==='all'||card.dataset.category===category)&&(status==='all'||card.dataset.review===status);card.hidden=!show;if(show)visible++;}}
 document.querySelector('#visible').textContent=`${{visible}} / ${{cards.length}} targets`;
}}
for(const id of ['search','category','status'])document.querySelector('#'+id).addEventListener('input',filter);
document.querySelector('#export').addEventListener('click',()=>{{
 const payload=cards.map(card=>({{target_id:card.dataset.targetId,...(saved[card.dataset.targetId]||{{status:'unreviewed',note:''}})}}));
 const blob=new Blob([JSON.stringify(payload,null,2)],{{type:'application/json'}});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='multistep-panel-review.json';a.click();URL.revokeObjectURL(a.href);
}});filter();
"""


def render_multistep_panel_html(
    cases: Iterable[MultistepPanelCase],
    *,
    title: str = "Familiar-target multistep retrosynthesis review",
    top_k: int = 3,
    metadata: Optional[Mapping[str, Any]] = None,
) -> str:
    """Render baseline and policy results as a self-contained HTML review."""

    values = tuple(cases)
    if not values:
        raise ValueError("multistep panel requires at least one case")
    if top_k < 1:
        raise ValueError("multistep panel top_k must be positive")
    solved_baseline = sum(bool(case.baseline.routes) for case in values)
    solved_policy = sum(bool(case.policy and case.policy.routes) for case in values)
    reference_count = sum(bool(case.reference_reactions) for case in values)
    partial_only_baseline = sum(
        not case.baseline.routes and bool(case.baseline.partial_routes)
        for case in values
    )
    policy_cases = tuple(case for case in values if case.policy is not None)
    changed = sum(
        _route_projection(case.baseline) != _route_projection(case.policy)
        for case in policy_cases
    )
    active = sum(bool(case.policy and case.policy.route_action_policy_active) for case in values)
    categories = sorted({case.target.category for case in values})
    category_options = "".join(
        f'<option value="{html.escape(value)}">{html.escape(value.title())}</option>'
        for value in categories
    )
    cards = "".join(
        _case_card(case, index, top_k=top_k)
        for index, case in enumerate(values, 1)
    )
    report_metadata = dict(metadata or {})
    storage_key = str(
        report_metadata.get("panel_id") or "multistep-familiar-panel-v1"
    )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body>'
        f'<header><div><h1>{html.escape(title)}</h1><p>Baseline and optional route-policy '
        'searches use identical chemistry validation and bounded-search settings. Review '
        'feasibility and strategic quality independently of structural validity.</p></div>'
        '<section class="summary">'
        f'<div><span>Targets</span><strong>{len(values)}</strong></div>'
        f'<div><span>Observed references</span><strong>{reference_count}</strong></div>'
        f'<div><span>Planner solved</span><strong>{solved_baseline}</strong></div>'
        f'<div><span>Planner partial-only</span><strong>{partial_only_baseline}</strong></div>'
        + (
            f'<div><span>Policy solved</span><strong>{solved_policy}</strong></div>'
            f'<div><span>Changed rankings</span><strong>{changed}</strong></div>'
            f'<div><span>Active policies</span><strong>{active}</strong></div>'
            if policy_cases
            else ""
        )
        + "</section></header>"
        '<main><section class="toolbar"><label>Search<input id="search" type="search" '
        'placeholder="Target, class, or SMILES"></label><label>Category<select id="category">'
        f'<option value="all">All</option>{category_options}</select></label>'
        '<label>Review status<select id="status"><option value="all">All</option>'
        '<option value="unreviewed">Unreviewed</option><option value="accept">Accept</option>'
        '<option value="question">Questionable</option><option value="reject">Reject</option>'
        '</select></label><button id="export" type="button">Export review JSON</button>'
        '<span class="visible" id="visible"></span></section>'
        f'<section>{cards}</section></main><script>{_script(storage_key)}</script></body></html>'
    )


def write_multistep_panel_artifacts(
    cases: Iterable[MultistepPanelCase],
    output_html: str | Path,
    *,
    output_json: Optional[str | Path] = None,
    title: str = "Familiar-target multistep retrosynthesis review",
    top_k: int = 3,
    metadata: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Write deterministic panel JSON and a self-contained review HTML."""

    values = tuple(cases)
    document = render_multistep_panel_html(
        values,
        title=title,
        top_k=top_k,
        metadata=metadata,
    )
    html_path = Path(output_html)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(document, encoding="utf-8", newline="\n")
    json_path = Path(output_json) if output_json is not None else None
    payload = {
        "schema_version": "1.0",
        "metadata": dict(metadata or {}),
        "cases": [case.to_dict() for case in values],
    }
    if json_path is not None:
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True),
            encoding="utf-8",
            newline="\n",
        )
    return {
        "target_count": len(values),
        "baseline_solved_count": sum(bool(case.baseline.routes) for case in values),
        "policy_solved_count": sum(bool(case.policy and case.policy.routes) for case in values),
        "changed_ranking_count": sum(
            _route_projection(case.baseline) != _route_projection(case.policy)
            for case in values
            if case.policy is not None
        ),
        "output_html": str(html_path.resolve()),
        "output_json": str(json_path.resolve()) if json_path is not None else None,
        "html_bytes": html_path.stat().st_size,
    }


__all__ = [
    "MultistepPanelCase",
    "MultistepPanelTarget",
    "load_multistep_panel_targets",
    "render_multistep_panel_html",
    "write_multistep_panel_artifacts",
]
