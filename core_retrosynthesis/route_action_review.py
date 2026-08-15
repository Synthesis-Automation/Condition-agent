"""Self-contained chemist review for promoted observed route-action labels."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import html
import json
from pathlib import Path
import random
from typing import Any, Iterable, Optional

from .html_report import reaction_svg
from .route_action_conversion import iter_route_action_evaluations
from .route_action_evaluation import (
    RouteActionEvaluation,
    RouteActionStepEvaluation,
)


DEFAULT_ROUTE_ACTION_REVIEW_SEED = 20_260_814

RouteActionReviewItem = tuple[RouteActionEvaluation, RouteActionStepEvaluation]


def is_promoted_route_action(step: RouteActionStepEvaluation) -> bool:
    """Return whether route evidence is usable beyond strict operator admission."""

    observed = step.observed_action
    return bool(observed.retained_edits_verified and not observed.operator_roundtrip_verified)


def _stratum(item: RouteActionReviewItem) -> tuple[str, str, str, str]:
    _, step = item
    observed = step.observed_action
    departing = (
        observed.departing_edit_descriptors[0]
        if observed.departing_edit_descriptors
        else "no_departing_edit"
    )
    return (
        observed.core_quality_status,
        observed.named_annotation or "unannotated",
        departing,
        "strategy" if observed.strategy_verified else "retained_only",
    )


def sample_promoted_route_actions(
    evaluations: Iterable[RouteActionEvaluation],
    *,
    sample_size: int = 100,
    seed: int = DEFAULT_ROUTE_ACTION_REVIEW_SEED,
) -> tuple[RouteActionReviewItem, ...]:
    """Return a deterministic round-robin sample across chemistry strata."""

    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    groups: dict[tuple[str, str, str, str], list[RouteActionReviewItem]] = defaultdict(list)
    for evaluation in evaluations:
        for step in evaluation.steps:
            item = (evaluation, step)
            if is_promoted_route_action(step):
                groups[_stratum(item)].append(item)
    if not groups:
        raise ValueError("no promoted route-action labels are available")
    for key, values in groups.items():
        values.sort(key=lambda item: item[1].evaluation_id)
        group_seed = int.from_bytes(
            hashlib.sha256(f"{seed}:{key}".encode("utf-8")).digest()[:8],
            byteorder="big",
        )
        random.Random(group_seed).shuffle(values)
    selected: list[RouteActionReviewItem] = []
    ordered_keys = sorted(groups)
    while len(selected) < sample_size:
        progressed = False
        for key in ordered_keys:
            values = groups[key]
            if values:
                selected.append(values.pop())
                progressed = True
                if len(selected) == sample_size:
                    break
        if not progressed:
            break
    return tuple(selected)


def _short(value: Optional[str], length: int = 32) -> str:
    if not value:
        return "unavailable"
    return value if len(value) <= length else f"{value[:length]}…"


def _badge(label: str, status: bool) -> str:
    state = "yes" if status else "no"
    return f'<span class="badge {state}">{html.escape(label)}</span>'


def _card(item: RouteActionReviewItem, index: int) -> str:
    route, step = item
    observed = step.observed_action
    route_id = route.source_route_id or route.tree_id
    descriptors = " · ".join(observed.departing_edit_descriptors) or "none"
    limitations = " · ".join(observed.limitations) or "none"
    annotation = observed.named_annotation or "unannotated"
    search = " ".join(
        (
            route_id,
            route.patent_id or "",
            annotation,
            descriptors,
            limitations,
            observed.target_smiles or "",
        )
    ).lower()
    return (
        f'<article class="action-card" data-search="{html.escape(search)}" '
        f'data-quality="{html.escape(observed.core_quality_status)}" '
        f'data-label="{"strategy" if observed.strategy_verified else "retained"}">'
        '<div class="card-heading"><div>'
        f'<strong>{html.escape(route_id)}</strong>'
        f'<small>{html.escape(route.patent_id or "unknown patent")} · depth '
        f'{step.retrosynthetic_depth} · {html.escape(step.step_id)}</small></div>'
        f'<span class="quality">{html.escape(observed.core_quality_status)}</span></div>'
        f'<div class="reaction">{reaction_svg(step.reaction_smiles)}</div>'
        '<div class="badges">'
        f'{_badge("product site", observed.product_site_verified)}'
        f'{_badge("retained edits", observed.retained_edits_verified)}'
        f'{_badge("synthon", observed.synthon_partition_verified)}'
        f'{_badge("exact precursors", observed.exact_precursors_verified)}'
        f'{_badge("strategy", observed.strategy_verified)}'
        f'{_badge("strict operator", observed.operator_roundtrip_verified)}'
        '</div>'
        '<dl>'
        f'<div><dt>Annotation</dt><dd>{html.escape(annotation)}</dd></div>'
        f'<div><dt>Retained OP1</dt><dd title="{html.escape(observed.retained_operator_id or "")}">'
        f'{html.escape(_short(observed.retained_operator_id))}</dd></div>'
        f'<div><dt>SITE1</dt><dd title="{html.escape(observed.disconnection_site_key or "")}">'
        f'{html.escape(_short(observed.disconnection_site_key))}</dd></div>'
        f'<div><dt>SYN1</dt><dd title="{html.escape(observed.synthon_signature or "")}">'
        f'{html.escape(_short(observed.synthon_signature))}</dd></div>'
        f'<div><dt>STRAT1</dt><dd title="{html.escape(observed.strategy_id or "")}">'
        f'{html.escape(_short(observed.strategy_id))}</dd></div>'
        f'<div><dt>Departing edits</dt><dd>{html.escape(descriptors)}</dd></div>'
        f'<div><dt>Strict rejection</dt><dd>{html.escape(observed.operator_admission_stage)}: '
        f'{html.escape(observed.operator_admission_reason or "none")}</dd></div>'
        f'<div><dt>Limitations</dt><dd>{html.escape(limitations)}</dd></div>'
        '</dl><div class="review"><label>Status<select class="review-status" '
        f'data-id="{html.escape(step.evaluation_id)}"><option value="unreviewed">Unreviewed</option>'
        '<option value="accept">Accept label</option><option value="reject">Reject label</option>'
        '<option value="uncertain">Uncertain</option></select></label>'
        '<label>Note<textarea class="review-note" '
        f'data-id="{html.escape(step.evaluation_id)}" placeholder="Chemistry issue or rationale"></textarea>'
        f'</label></div><input type="hidden" value="{index}"></article>'
    )


def _stylesheet() -> str:
    return """
:root{--bg:#f4f1ea;--paper:#fffdf8;--ink:#20231f;--muted:#676b63;--line:#d8d2c5;
--good:#1e6b4d;--bad:#9c3f35;--accent:#245f73}*{box-sizing:border-box}body{margin:0;background:var(--bg);
color:var(--ink);font:14px/1.45 system-ui,sans-serif}header,main{max-width:1500px;margin:auto;padding:22px}
header p{color:var(--muted);max-width:850px}.summary,.toolbar{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:16px}
.summary>div{background:var(--paper);border:1px solid var(--line);border-radius:8px;padding:10px 16px}
.summary span{display:block;color:var(--muted);font-size:11px;text-transform:uppercase}.summary strong{font-size:20px}
.toolbar{position:sticky;top:0;z-index:3;background:rgba(244,241,234,.96);padding:10px 0}.toolbar label{display:flex;
align-items:center;gap:6px}.toolbar input,.toolbar select,button,textarea,.review select{border:1px solid var(--line);
border-radius:6px;background:white;padding:7px}button{cursor:pointer}.action-card{background:var(--paper);border:1px solid var(--line);
border-radius:10px;padding:14px;margin:0 0 18px}.card-heading{display:flex;justify-content:space-between;gap:10px}
.card-heading small{display:block;color:var(--muted)}.quality{border:1px solid var(--line);border-radius:20px;padding:4px 10px;
height:max-content}.reaction{overflow-x:auto;background:white;margin:12px 0;border-radius:6px}.reaction svg{min-width:700px;width:100%;height:250px}
.badges{display:flex;gap:7px;flex-wrap:wrap}.badge{padding:3px 8px;border-radius:12px;font-size:11px}.badge.yes{background:#dceee5;
color:var(--good)}.badge.no{background:#f0dedb;color:var(--bad)}dl{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:4px 18px}
dl div{border-top:1px solid var(--line);padding:6px 0;min-width:0}dt{color:var(--muted);font-size:11px}dd{margin:2px 0;
overflow-wrap:anywhere}.review{display:grid;grid-template-columns:220px 1fr;gap:12px;background:#f7f4ed;padding:10px;border-radius:7px}
.review label{display:flex;flex-direction:column;gap:4px}.review textarea{min-height:64px}.hidden{display:none}@media(max-width:800px){
dl{grid-template-columns:1fr}.review{grid-template-columns:1fr}.reaction svg{min-width:560px}}
"""


def _script(report_key: str, metadata: dict[str, Any]) -> str:
    return f"""
const key={json.dumps(report_key)};const metadata={json.dumps(metadata, sort_keys=True)};
const load=()=>JSON.parse(localStorage.getItem(key)||'{{}}');let reviews=load();
document.querySelectorAll('.review-status').forEach(el=>{{el.value=(reviews[el.dataset.id]||{{}}).status||'unreviewed';
el.onchange=save;}});document.querySelectorAll('.review-note').forEach(el=>{{el.value=(reviews[el.dataset.id]||{{}}).note||'';el.oninput=save;}});
function save(e){{const id=e.target.dataset.id;reviews[id]=reviews[id]||{{}};if(e.target.classList.contains('review-status'))reviews[id].status=e.target.value;
else reviews[id].note=e.target.value;localStorage.setItem(key,JSON.stringify(reviews));}}
function filter(){{const q=document.querySelector('#search').value.toLowerCase();const quality=document.querySelector('#quality').value;
const label=document.querySelector('#label').value;let visible=0;document.querySelectorAll('.action-card').forEach(card=>{{const show=card.dataset.search.includes(q)&&
(quality==='all'||card.dataset.quality===quality)&&(label==='all'||card.dataset.label===label);card.classList.toggle('hidden',!show);visible+=show?1:0;}});
document.querySelector('#visible').textContent=visible+' visible';}}document.querySelectorAll('.toolbar input,.toolbar select').forEach(el=>el.oninput=filter);filter();
document.querySelector('#export').onclick=()=>{{const payload={{metadata,reviews}};const blob=new Blob([JSON.stringify(payload,null,2)],{{type:'application/json'}});
const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='route-action-review.json';a.click();URL.revokeObjectURL(a.href);}};
"""


def render_route_action_review_html(
    evaluations: Iterable[RouteActionEvaluation],
    *,
    sample_size: int = 100,
    seed: int = DEFAULT_ROUTE_ACTION_REVIEW_SEED,
    title: str = "Promoted route-action label review",
) -> str:
    """Render a stratified, self-contained review page."""

    materialized = tuple(evaluations)
    promoted_count = sum(
        is_promoted_route_action(step)
        for route in materialized
        for step in route.steps
    )
    sampled = sample_promoted_route_actions(
        materialized,
        sample_size=min(sample_size, promoted_count),
        seed=seed,
    )
    strategy_count = sum(item[1].observed_action.strategy_verified for item in sampled)
    departing_count = sum(bool(item[1].observed_action.departing_edit_count) for item in sampled)
    cards = "".join(_card(item, index) for index, item in enumerate(sampled, 1))
    metadata = {
        "seed": seed,
        "sample_size": len(sampled),
        "evaluation_ids": [item[1].evaluation_id for item in sampled],
    }
    report_key = f"route-action-review:{seed}:{hashlib.sha256('|'.join(metadata['evaluation_ids']).encode()).hexdigest()[:16]}"
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body><header>'
        f'<h1>{html.escape(title)}</h1><p>Observed one-step route evidence promoted for learning even '
        'when the stricter executable operator cannot be admitted. Review the chemistry label; acceptance '
        'here does not admit an executable template.</p></header><main><section class="summary">'
        f'<div><span>Promoted population</span><strong>{promoted_count}</strong></div>'
        f'<div><span>Reviewed sample</span><strong>{len(sampled)}</strong></div>'
        f'<div><span>STRAT1 labels</span><strong>{strategy_count}</strong></div>'
        f'<div><span>Departing evidence</span><strong>{departing_count}</strong></div></section>'
        '<section class="toolbar"><label>Search<input id="search" type="search" placeholder="Route, patent, chemistry"></label>'
        '<label>Core quality<select id="quality"><option value="all">All</option><option value="pass">Pass</option>'
        '<option value="review">Review</option></select></label><label>Label<select id="label"><option value="all">All</option>'
        '<option value="strategy">STRAT1</option><option value="retained">Retained-only</option></select></label>'
        '<button id="export" type="button">Export review</button><span id="visible"></span></section>'
        f'<section>{cards}</section></main><script>{_script(report_key, metadata)}</script></body></html>'
    )


def write_route_action_review_html(
    source_evaluations: str | Path,
    output_html: str | Path,
    *,
    sample_size: int = 100,
    seed: int = DEFAULT_ROUTE_ACTION_REVIEW_SEED,
    title: str = "Promoted route-action label review",
) -> dict[str, Any]:
    """Load route-action labels and write a self-contained chemist review."""

    evaluations = tuple(iter_route_action_evaluations(source_evaluations))
    promoted = [
        (route, step)
        for route in evaluations
        for step in route.steps
        if is_promoted_route_action(step)
    ]
    document = render_route_action_review_html(
        evaluations,
        sample_size=sample_size,
        seed=seed,
        title=title,
    )
    output = Path(output_html)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8", newline="\n")
    sampled = sample_promoted_route_actions(
        evaluations,
        sample_size=min(sample_size, len(promoted)),
        seed=seed,
    )
    return {
        "source_evaluations": str(Path(source_evaluations).resolve()),
        "output_html": str(output.resolve()),
        "promoted_step_count": len(promoted),
        "sample_size": len(sampled),
        "seed": seed,
        "evaluation_ids": [item[1].evaluation_id for item in sampled],
        "html_bytes": output.stat().st_size,
    }


__all__ = [
    "DEFAULT_ROUTE_ACTION_REVIEW_SEED",
    "is_promoted_route_action",
    "render_route_action_review_html",
    "sample_promoted_route_actions",
    "write_route_action_review_html",
]
