"""Self-contained chemist review for minimized route-core projections."""

from __future__ import annotations

import html
import json
from pathlib import Path
import random
from typing import Any, Iterable, Mapping, Optional

from .html_report import molecule_svg, reaction_svg
from .route_core import RouteCoreLineageLink, RouteCoreProjection, RouteCoreStep
from .route_core_conversion import iter_route_core_projections


DEFAULT_ROUTE_CORE_REVIEW_SEED = 20_260_814


def sample_route_core_projections(
    projections: Iterable[RouteCoreProjection],
    *,
    sample_size: int = 50,
    seed: int = DEFAULT_ROUTE_CORE_REVIEW_SEED,
) -> tuple[RouteCoreProjection, ...]:
    """Return a deterministic simple-random route-core sample."""

    ordered = sorted(
        projections,
        key=lambda item: item.source_route_id or item.route_core_id,
    )
    if sample_size < 1 or sample_size > len(ordered):
        raise ValueError(f"sample_size must be between 1 and {len(ordered)}")
    return tuple(random.Random(seed).sample(ordered, sample_size))


def _short(value: Optional[str], length: int = 18) -> str:
    if not value:
        return "unavailable"
    return value if len(value) <= length else f"{value[:length]}…"


def _lineage_connector(
    link: Optional[RouteCoreLineageLink],
) -> str:
    if link is None:
        status = "branch"
        label = "convergent branch"
        detail = "events are not a single direct lineage"
    else:
        status = link.status
        label = {
            "unique": "unique atom lineage",
            "symmetry_ambiguous": "symmetry-equivalent lineage",
            "bounded_ambiguous": "bounded lineage ambiguity",
        }.get(link.status, "unresolved lineage")
        detail = f"{link.candidate_count} correspondence candidate(s)"
    return (
        f'<div class="core-connector {html.escape(status)}">'
        '<div class="core-arrow"><span class="forward-arrow">⟶</span>'
        '<span class="retro-arrow">⟵</span></div>'
        f'<strong>{html.escape(label)}</strong>'
        f'<small>{html.escape(detail)}</small></div>'
    )


def _step_card(step: RouteCoreStep, index: int, total: int) -> str:
    reaction = step.render_reaction_smiles or step.minimum_reaction_smiles or ""
    title = "Target-forming core event" if index == total else f"Core event {index}"
    edits = " · ".join(step.edit_tokens) or "no resolved edit tokens"
    transformation = step.transformation_class or "unresolved transformation"
    return (
        f'<article class="core-step {html.escape(step.quality_status)}">'
        '<div class="core-step-heading">'
        f'<strong>{html.escape(title)}</strong>'
        f'<span class="quality-tag">{html.escape(step.quality_status)}</span></div>'
        f'<div class="core-reaction">{reaction_svg(reaction)}</div>'
        f'<div class="core-class">{html.escape(transformation)}</div>'
        f'<div class="edit-tokens">{html.escape(edits)}</div>'
        '<dl><div><dt>Active atoms</dt>'
        f'<dd>{step.active_atom_count}</dd></div>'
        '<div><dt>Events</dt>'
        f'<dd>{step.event_count}</dd></div>'
        '<div><dt>R groups</dt>'
        f'<dd>{step.display_r_group_count}</dd></div></dl>'
        f'<code title="{html.escape(step.reaction_core_id or "")}">'
        f'{html.escape(_short(step.reaction_core_id))}</code>'
        f'<code>{html.escape(step.minimum_reaction_smiles or "")}</code>'
        '</article>'
    )


def _route_card(projection: RouteCoreProjection, index: int) -> str:
    route_id = projection.source_route_id or projection.route_core_id
    link_by_pair = {
        (link.producer_reaction_node_id, link.consumer_reaction_node_id): link
        for link in projection.lineage_links
    }
    track: list[str] = []
    for position, step in enumerate(projection.steps, 1):
        track.append(_step_card(step, position, len(projection.steps)))
        if position < len(projection.steps):
            following = projection.steps[position]
            track.append(
                _lineage_connector(
                    link_by_pair.get(
                        (step.reaction_node_id, following.reaction_node_id)
                    )
                )
            )
    quality = "resolved" if projection.fully_chemistry_resolved else "partial"
    lineage = "connected" if projection.fully_lineage_connected else "partial"
    search = " ".join(
        (
            route_id,
            projection.patent_id or "",
            projection.split or "",
            projection.target_smiles,
            quality,
            lineage,
        )
    ).lower()
    return (
        f'<details class="route-core-card" data-route-id="{html.escape(route_id)}" '
        f'data-quality="{quality}" data-lineage="{lineage}" '
        f'data-search="{html.escape(search)}"{ " open" if index <= 2 else ""}>'
        '<summary>'
        f'<span class="route-index">{index}</span>'
        '<span class="route-title">'
        f'<strong>{html.escape(route_id)}</strong>'
        f'<small>{html.escape(projection.patent_id or "unknown patent")}</small>'
        '</span>'
        f'<span class="tag">{projection.reaction_count} core events</span>'
        f'<span class="tag {quality}">{quality} chemistry</span>'
        f'<span class="tag {lineage}">{lineage} lineage</span>'
        '</summary><div class="route-core-body">'
        '<section class="core-route-heading"><div>'
        '<h2>Minimized route core</h2>'
        '<p><span class="forward-copy">Forward synthesis order</span>'
        '<span class="retro-copy">Retrosynthetic order</span></p></div>'
        '<div class="target-anchor"><span>Full target anchor</span>'
        f'{molecule_svg(projection.target_smiles, width=280, height=170)}</div>'
        '</section>'
        f'<div class="core-scroll"><div class="core-track">{"".join(track)}</div></div>'
        '<section class="route-core-identity">'
        f'<div><span>Exact route core</span><code>{html.escape(projection.exact_route_core_key)}</code></div>'
        f'<div><span>Typed route core</span><code>{html.escape(projection.typed_route_core_key)}</code></div>'
        f'<div><span>Shape route core</span><code>{html.escape(projection.shape_route_core_key)}</code></div>'
        '</section>'
        '<section class="review-fields">'
        f'<label for="status-{index}">Review</label>'
        f'<select id="status-{index}" class="review-status" data-route="{html.escape(route_id)}">'
        '<option value="unreviewed">Unreviewed</option>'
        '<option value="accept">Accept</option>'
        '<option value="question">Questionable</option>'
        '<option value="reject">Reject</option></select>'
        f'<label for="note-{index}">Notes</label>'
        f'<textarea id="note-{index}" class="review-note" data-route="{html.escape(route_id)}" '
        'rows="2" placeholder="Core, context, or lineage note"></textarea>'
        '</section></div></details>'
    )


def _stylesheet() -> str:
    return """
:root { color-scheme:light; --ink:#17242c; --muted:#68777e; --line:#d7e0e3;
  --panel:#fff; --page:#edf3f4; --accent:#17677a; --soft:#f7fafb;
  --good:#287552; --warn:#a16b10; --bad:#a33a3a; }
* { box-sizing:border-box; }
body { margin:0; color:var(--ink); background:var(--page);
  font:14px/1.42 Inter,Segoe UI,Arial,sans-serif; }
header { padding:24px max(22px,calc((100vw - 1600px)/2)); color:#fff;
  background:linear-gradient(120deg,#143b49,#17677a); }
header h1 { margin:0 0 4px; font-size:27px; font-weight:520; }
header p { margin:0; opacity:.88; }
main { max-width:1600px; margin:auto; padding:18px; }
.summary { display:grid; grid-template-columns:repeat(5,1fr); gap:9px; margin-bottom:12px; }
.summary div { padding:10px 13px; background:#fff; border:1px solid var(--line); border-radius:8px; }
.summary span { color:var(--muted); font-size:11px; }.summary strong { display:block; font-size:19px; }
.toolbar { position:sticky; top:0; z-index:5; display:flex; flex-wrap:wrap; gap:9px;
  align-items:end; margin-bottom:12px; padding:10px 12px; background:#fff;
  border:1px solid var(--line); border-radius:8px; }
.toolbar label { display:grid; gap:3px; color:var(--muted); font-size:11px; }
.toolbar input,.toolbar select,.toolbar button,.review-fields select,.review-fields textarea {
  padding:7px 9px; color:var(--ink); background:#fff; border:1px solid #b9c7cc;
  border-radius:6px; font:inherit; }.toolbar button { cursor:pointer; background:var(--soft); }
#visible-count { margin-left:auto; padding-bottom:7px; color:var(--muted); }
.route-core-card { margin-bottom:11px; overflow:hidden; background:#fff;
  border:1px solid var(--line); border-radius:9px; }
.route-core-card>summary { display:flex; gap:8px; align-items:center; padding:11px 13px;
  cursor:pointer; background:var(--soft); }.route-index { display:grid; place-items:center;
  width:27px; height:27px; color:#124d5c; background:#dceef2; border-radius:50%; }
.route-title { display:grid; margin-right:auto; }.route-title small { color:var(--muted); }
.tag { padding:3px 7px; color:#31535d; background:#e6eff2; border-radius:999px; font-size:11px; }
.tag.resolved,.tag.connected { color:#215f43; background:#e2f1e9; }
.tag.partial { color:#855b12; background:#f7edd7; }
.route-core-body { padding:14px; }.core-route-heading { display:flex; align-items:center;
  justify-content:space-between; gap:15px; }.core-route-heading h2 { margin:0; font-size:18px; }
.core-route-heading p { margin:3px 0 0; color:var(--muted); }.target-anchor { display:flex;
  align-items:center; gap:8px; color:var(--muted); font-size:11px; }.target-anchor svg { width:220px; height:130px; }
.core-scroll { overflow-x:auto; padding:10px 2px; }.core-track { display:flex;
  width:max-content; min-width:100%; align-items:stretch; }.core-step { display:grid;
  flex:0 0 440px; padding:9px; border:1px solid var(--line); border-top:4px solid var(--warn);
  border-radius:8px; background:#fff; }.core-step.pass { border-top-color:var(--good); }
.core-step.blocked,.core-step.unavailable { border-top-color:var(--bad); }
.core-step-heading { display:flex; justify-content:space-between; gap:8px; }
.quality-tag { color:var(--muted); font-size:11px; }.core-reaction svg { width:100%; height:220px; }
.core-class { font-weight:600; text-align:center; }.edit-tokens { min-height:38px; color:#43545b;
  font:11px/1.35 Consolas,monospace; text-align:center; overflow-wrap:anywhere; }
.core-step dl { display:grid; grid-template-columns:repeat(3,1fr); margin:6px 0; }
.core-step dl div { text-align:center; }.core-step dt { color:var(--muted); font-size:10px; }
.core-step dd { margin:0; }.core-step code { max-height:40px; overflow:auto; }
.core-connector { display:grid; flex:0 0 165px; align-content:center; justify-items:center;
  padding:8px; color:#31535d; text-align:center; }.core-connector strong { font-size:11px; }
.core-connector small { color:var(--muted); }.core-arrow { color:var(--accent); font:42px/1 Arial; }
.core-connector.symmetry_ambiguous { color:var(--warn); }.core-connector.branch,
.core-connector.component_unresolved,.core-connector.graph_mismatch { color:var(--bad); }
.retro-arrow,.retro-copy { display:none; }body.retro-direction .core-track { flex-direction:row-reverse; }
body.retro-direction .forward-arrow,body.retro-direction .forward-copy { display:none; }
body.retro-direction .retro-arrow,body.retro-direction .retro-copy { display:inline; }
.route-core-identity { display:grid; grid-template-columns:repeat(3,1fr); gap:7px; margin-top:5px; }
.route-core-identity div { min-width:0; padding:7px; background:var(--soft); border-radius:6px; }
.route-core-identity span { color:var(--muted); font-size:10px; }.route-core-identity code,
.core-step code { display:block; overflow-wrap:anywhere; color:#405057; font:10px/1.3 Consolas,monospace; }
.review-fields { display:grid; grid-template-columns:auto 180px auto 1fr; gap:7px;
  align-items:center; margin-top:10px; padding-top:10px; border-top:1px solid var(--line); }
.hidden { display:none; }.render-error { padding:50px; color:var(--bad); text-align:center; }
@media (max-width:900px) { .summary { grid-template-columns:1fr 1fr; }.toolbar { position:static; }
  .core-route-heading { align-items:start; flex-direction:column; }.target-anchor { display:none; }
  .route-core-identity { grid-template-columns:1fr; }.review-fields { grid-template-columns:1fr; } }
"""


def _script(report_key: str, metadata: Mapping[str, Any]) -> str:
    key = json.dumps(report_key)
    payload = json.dumps(metadata, sort_keys=True).replace("</", "<\\/")
    return f"""
const reportKey = {key}; const sampleMetadata = {payload};
const cards = [...document.querySelectorAll('.route-core-card')];
const search = document.getElementById('search-filter');
const quality = document.getElementById('quality-filter');
const lineage = document.getElementById('lineage-filter');
const visible = document.getElementById('visible-count');
function readReviews() {{ try {{ return JSON.parse(localStorage.getItem(reportKey) || '{{}}'); }} catch (_) {{ return {{}}; }} }}
let reviews = readReviews();
function filterCards() {{ const query = search.value.trim().toLowerCase(); let count = 0;
  cards.forEach(card => {{ const show = (!query || card.dataset.search.includes(query)) &&
    (quality.value === 'all' || card.dataset.quality === quality.value) &&
    (lineage.value === 'all' || card.dataset.lineage === lineage.value);
    card.classList.toggle('hidden', !show); if (show) count += 1; }});
  visible.textContent = `${{count}} / ${{cards.length}} routes`; }}
cards.forEach(card => {{ const id = card.dataset.routeId; const saved = reviews[id] || {{status:'unreviewed',note:''}};
  card.querySelector('.review-status').value = saved.status; card.querySelector('.review-note').value = saved.note;
  card.querySelectorAll('.review-status,.review-note').forEach(control => control.addEventListener('change', () => {{
    reviews[id] = {{status:card.querySelector('.review-status').value,note:card.querySelector('.review-note').value}};
    try {{ localStorage.setItem(reportKey, JSON.stringify(reviews)); }} catch (_) {{}} }})); }});
[search,quality,lineage].forEach(control => control.addEventListener('input', filterCards));
document.getElementById('toggle-direction').addEventListener('click', event => {{ const retro = document.body.classList.toggle('retro-direction');
  event.currentTarget.textContent = retro ? 'Show forward synthesis' : 'Show retrosynthesis'; }});
document.getElementById('expand-all').addEventListener('click', () => cards.forEach(card => card.open = true));
document.getElementById('collapse-all').addEventListener('click', () => cards.forEach(card => card.open = false));
document.getElementById('export-review').addEventListener('click', () => {{ const data = {{sample:sampleMetadata,reviews}};
  const blob = new Blob([JSON.stringify(data,null,2)],{{type:'application/json'}}); const link = document.createElement('a');
  link.href=URL.createObjectURL(blob); link.download='route-core-review.json'; link.click(); URL.revokeObjectURL(link.href); }});
filterCards();
"""


def render_route_core_review_html(
    projections: Iterable[RouteCoreProjection],
    *,
    sample_size: int = 50,
    seed: int = DEFAULT_ROUTE_CORE_REVIEW_SEED,
    title: str = "Minimized route-core review",
) -> str:
    """Render route-core projections as a self-contained chemistry review."""

    sampled = sample_route_core_projections(
        projections,
        sample_size=sample_size,
        seed=seed,
    )
    total_steps = sum(item.reaction_count for item in sampled)
    resolved = sum(item.fully_chemistry_resolved for item in sampled)
    connected = sum(item.fully_lineage_connected for item in sampled)
    motif_count = len(
        {key for item in sampled for key in item.typed_motif_keys}
    )
    cards = "".join(
        _route_card(item, index) for index, item in enumerate(sampled, 1)
    )
    metadata = {
        "seed": seed,
        "sample_size": sample_size,
        "route_core_ids": [item.route_core_id for item in sampled],
        "source_route_ids": [item.source_route_id for item in sampled],
    }
    report_key = f"route-core-review:{seed}:{sample_size}"
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body>'
        f'<header><h1>{html.escape(title)}</h1><p>Context-preserving minimized reactions; '
        'projection evidence, not executable multistep templates.</p></header><main>'
        '<section class="summary">'
        f'<div><span>Routes</span><strong>{sample_size}</strong></div>'
        f'<div><span>Core events</span><strong>{total_steps}</strong></div>'
        f'<div><span>Chemistry resolved</span><strong>{resolved}</strong></div>'
        f'<div><span>Lineage connected</span><strong>{connected}</strong></div>'
        f'<div><span>Typed motifs</span><strong>{motif_count}</strong></div></section>'
        '<section class="toolbar"><label>Search<input id="search-filter" type="search" '
        'placeholder="Route, patent, target"></label><label>Chemistry<select id="quality-filter">'
        '<option value="all">All</option><option value="resolved">Resolved</option>'
        '<option value="partial">Partial</option></select></label><label>Lineage<select id="lineage-filter">'
        '<option value="all">All</option><option value="connected">Connected</option>'
        '<option value="partial">Partial</option></select></label>'
        '<button id="toggle-direction" type="button">Show retrosynthesis</button>'
        '<button id="expand-all" type="button">Expand all</button>'
        '<button id="collapse-all" type="button">Collapse all</button>'
        '<button id="export-review" type="button">Export review</button>'
        '<span id="visible-count"></span></section>'
        f'<section id="routes">{cards}</section></main>'
        f'<script>{_script(report_key, metadata)}</script></body></html>'
    )


def write_route_core_review_html(
    source_route_cores: str | Path,
    output_html: str | Path,
    *,
    sample_size: int = 50,
    seed: int = DEFAULT_ROUTE_CORE_REVIEW_SEED,
    title: str = "Minimized route-core review",
) -> dict[str, Any]:
    """Load, sample and write a self-contained route-core review."""

    projections = tuple(iter_route_core_projections(source_route_cores))
    document = render_route_core_review_html(
        projections,
        sample_size=sample_size,
        seed=seed,
        title=title,
    )
    output = Path(output_html)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8", newline="\n")
    sampled = sample_route_core_projections(
        projections,
        sample_size=sample_size,
        seed=seed,
    )
    return {
        "source_route_cores": str(Path(source_route_cores).resolve()),
        "output_html": str(output.resolve()),
        "sample_size": sample_size,
        "seed": seed,
        "route_core_ids": [item.route_core_id for item in sampled],
        "source_route_ids": [item.source_route_id for item in sampled],
        "html_bytes": output.stat().st_size,
    }


__all__ = [
    "DEFAULT_ROUTE_CORE_REVIEW_SEED",
    "render_route_core_review_html",
    "sample_route_core_projections",
    "write_route_core_review_html",
]
