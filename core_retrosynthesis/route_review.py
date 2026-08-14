"""Self-contained visual review for curated retrosynthesis routes."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import random
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

from .html_report import molecule_svg, reaction_svg


DEFAULT_ROUTE_REVIEW_SEED = 20_260_814


def load_curated_routes(path: str | Path) -> tuple[dict[str, Any], ...]:
    """Load curated route JSONL from a plain or gzip-compressed file."""

    source = Path(path)
    records: list[dict[str, Any]] = []
    if source.suffix.lower() == ".gz":
        handle = gzip.open(source, "rt", encoding="utf-8")
    else:
        handle = source.open("r", encoding="utf-8")
    with handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"Route record {line_number} is not an object")
            records.append(value)
    return tuple(records)


def sample_route_records(
    records: Iterable[Mapping[str, Any]],
    *,
    sample_size: int = 50,
    seed: int = DEFAULT_ROUTE_REVIEW_SEED,
) -> tuple[Mapping[str, Any], ...]:
    """Return a reproducible simple-random sample independent of file order."""

    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    ordered = sorted(records, key=lambda row: str(row.get("route_id") or ""))
    if sample_size > len(ordered):
        raise ValueError(
            f"sample_size {sample_size} exceeds {len(ordered)} route records"
        )
    return tuple(random.Random(seed).sample(ordered, sample_size))


def _route_card(route: Mapping[str, Any], index: int) -> str:
    route_id = str(route.get("route_id") or f"route-{index}")
    patent_id = str(route.get("patent_id") or "unknown")
    split = str(route.get("split") or "unknown")
    target = str(route.get("target_smiles") or "")
    steps = tuple(route.get("steps") or ())
    reduction = int(route.get("abstraction_reduction") or 0)
    step_blocks = []
    for fallback_position, step in enumerate(steps):
        position = int(step.get("retrosynthetic_position", fallback_position))
        original = str(step.get("reaction_smiles") or "")
        abstracted = str(step.get("abstracted_reaction_smiles") or "")
        reagents = str(step.get("reagents_smiles") or "")
        source_id = str(step.get("source_reaction_id") or "unknown")
        internal = ", ".join(step.get("internal_precursor_smiles") or ()) or "none"
        terminal = ", ".join(step.get("terminal_precursor_smiles") or ()) or "none"
        abstraction = (
            '<div class="reaction-panel abstracted">'
            '<h4>Higher-level abstraction</h4>'
            f'<div class="reaction-svg">{reaction_svg(abstracted)}</div>'
            f'<code>{html.escape(abstracted)}</code></div>'
            if abstracted
            else (
                '<div class="reaction-panel abstracted absorbed">'
                '<h4>Higher-level abstraction</h4>'
                '<p>Absorbed into an adjacent higher-level transformation.</p>'
                "</div>"
            )
        )
        step_blocks.append(
            '<section class="route-step">'
            '<div class="step-heading">'
            f'<span class="step-number">{position + 1}</span>'
            f'<h3>Retrosynthetic step {position + 1}</h3>'
            f'<span class="source-id">source {html.escape(source_id)}</span>'
            "</div>"
            '<div class="reaction-comparison">'
            '<div class="reaction-panel original">'
            '<h4>Original recorded reaction</h4>'
            f'<div class="reaction-svg">{reaction_svg(original)}</div>'
            f'<code>{html.escape(original)}</code></div>'
            f"{abstraction}</div>"
            '<dl class="step-evidence">'
            f'<div><dt>Reagents</dt><dd>{html.escape(reagents or "none recorded")}</dd></div>'
            f'<div><dt>Internal precursor</dt><dd>{html.escape(internal)}</dd></div>'
            f'<div><dt>Terminal precursors</dt><dd>{html.escape(terminal)}</dd></div>'
            "</dl></section>"
        )
    search_value = " ".join((route_id, patent_id, split, target)).lower()
    return (
        f'<details class="route-card" data-route-id="{html.escape(route_id)}" '
        f'data-split="{html.escape(split)}" data-steps="{len(steps)}" '
        f'data-status="unreviewed" data-search="{html.escape(search_value)}"'
        f"{' open' if index <= 2 else ''}>"
        "<summary>"
        f'<span class="route-index">{index}</span>'
        '<span class="route-title">'
        f'<strong>{html.escape(route_id)}</strong>'
        f'<span>{html.escape(patent_id)}</span></span>'
        f'<span class="tag">{html.escape(split)}</span>'
        f'<span class="tag">{len(steps)} steps</span>'
        f'<span class="tag">reduction {reduction}</span>'
        "</summary>"
        '<div class="route-body">'
        '<section class="target-panel">'
        '<div><h2>Target</h2>'
        f'<code>{html.escape(target)}</code></div>'
        f'<div class="target-svg">{molecule_svg(target, width=420, height=240)}</div>'
        '<div class="review-fields">'
        f'<label for="status-{index}">Review status</label>'
        f'<select id="status-{index}" class="review-status" data-route="{html.escape(route_id)}">'
        '<option value="unreviewed">Unreviewed</option>'
        '<option value="accept">Accept</option>'
        '<option value="question">Questionable</option>'
        '<option value="reject">Reject</option></select>'
        f'<label for="note-{index}">Notes</label>'
        f'<textarea id="note-{index}" class="review-note" data-route="{html.escape(route_id)}" '
        'rows="3" placeholder="Optional chemistry or data-quality note"></textarea>'
        "</div></section>"
        f"{''.join(step_blocks)}</div></details>"
    )


def _stylesheet() -> str:
    return """
:root { color-scheme:light; --ink:#17242c; --muted:#64727a; --line:#d8e0e4;
  --panel:#ffffff; --page:#eef3f5; --soft:#f7f9fa; --accent:#17677a;
  --accent-soft:#dceef2; --good:#287552; --warn:#9a6813; --bad:#a33a3a; }
* { box-sizing:border-box; }
body { margin:0; color:var(--ink); background:var(--page);
  font:14px/1.45 Inter,Segoe UI,Arial,sans-serif; }
header { padding:26px max(24px,calc((100vw - 1500px)/2)); color:#fff;
  background:linear-gradient(120deg,#143b49,#17677a); }
header h1 { margin:0 0 5px; font-size:27px; font-weight:500; }
header p { margin:0; opacity:.88; }
main { max-width:1500px; margin:auto; padding:20px; }
.summary { display:grid; grid-template-columns:repeat(4,minmax(120px,1fr));
  gap:10px; margin-bottom:14px; }
.summary div { padding:12px 14px; background:var(--panel); border:1px solid var(--line);
  border-radius:9px; }
.summary span { display:block; color:var(--muted); font-size:12px; }
.summary strong { display:block; margin-top:2px; font-size:20px; font-weight:500; }
.toolbar { position:sticky; top:0; z-index:5; display:flex; flex-wrap:wrap; gap:10px;
  align-items:end; padding:11px 13px; margin-bottom:14px; background:var(--panel);
  border:1px solid var(--line); border-radius:9px; box-shadow:0 2px 9px #18303912; }
.toolbar label { display:grid; gap:3px; color:var(--muted); font-size:12px; }
.toolbar input,.toolbar select,.review-fields select,.review-fields textarea {
  padding:7px 9px; color:var(--ink); background:#fff; border:1px solid #b8c6cc;
  border-radius:6px; font:inherit; }
.toolbar button { padding:8px 11px; color:var(--ink); background:var(--soft);
  border:1px solid #b8c6cc; border-radius:6px; cursor:pointer; }
#visible-count { margin-left:auto; padding-bottom:7px; color:var(--muted); }
.route-card { margin-bottom:12px; overflow:hidden; background:var(--panel);
  border:1px solid var(--line); border-radius:10px; }
.route-card[data-status="accept"] { border-left:5px solid var(--good); }
.route-card[data-status="question"] { border-left:5px solid var(--warn); }
.route-card[data-status="reject"] { border-left:5px solid var(--bad); }
.route-card>summary { display:flex; flex-wrap:wrap; gap:9px; align-items:center;
  padding:12px 14px; cursor:pointer; background:var(--soft); }
.route-index,.step-number { display:grid; place-items:center; flex:0 0 auto; width:27px;
  height:27px; color:#124d5c; background:var(--accent-soft); border-radius:50%; font-weight:500; }
.route-title { display:grid; margin-right:auto; }.route-title span { color:var(--muted); font-size:12px; }
.tag,.source-id { padding:3px 7px; color:#31535d; background:#e6eff2;
  border-radius:999px; font-size:11px; }
.route-body { padding:15px; }.target-panel { display:grid;
  grid-template-columns:minmax(220px,.6fr) minmax(320px,1fr) minmax(220px,.55fr);
  gap:14px; align-items:center; padding-bottom:15px; }
.target-panel h2,.route-step h3,.reaction-panel h4 { margin:0 0 7px; font-weight:500; }
.target-svg svg,.reaction-svg svg { display:block; width:100%; height:auto; max-height:250px; }
.review-fields { display:grid; gap:5px; }.review-fields textarea { resize:vertical; }
.route-step { padding:14px 0; border-top:1px solid var(--line); }
.step-heading { display:flex; gap:9px; align-items:center; margin-bottom:10px; }
.step-heading h3 { margin:0; }.source-id { margin-left:auto; }
.reaction-comparison { display:grid; grid-template-columns:1fr 1fr; gap:12px; }
.reaction-panel { min-width:0; padding:10px; background:var(--soft); border:1px solid var(--line);
  border-radius:8px; }.reaction-panel.absorbed { display:grid; place-items:center; min-height:190px;
  color:var(--muted); text-align:center; }.reaction-panel.absorbed h4 { align-self:start; justify-self:start; }
.step-evidence { display:grid; grid-template-columns:1fr 1fr 1fr; gap:8px; margin:9px 0 0; }
.step-evidence div { min-width:0; padding:7px 9px; background:#f5f8f9; border-radius:6px; }
.step-evidence dt { color:var(--muted); font-size:11px; }.step-evidence dd { margin:2px 0 0;
  overflow-wrap:anywhere; font-family:Consolas,monospace; font-size:12px; }
code { display:block; overflow-wrap:anywhere; white-space:normal; color:#35454d;
  font:11px/1.35 Consolas,monospace; }
.hidden { display:none; }.render-error { padding:35px; color:var(--bad); text-align:center; }
@media (max-width:900px) { .summary { grid-template-columns:1fr 1fr; }
  .target-panel,.reaction-comparison { grid-template-columns:1fr; }
  .step-evidence { grid-template-columns:1fr; }.toolbar { position:static; } }
"""


def _script(report_key: str, sample_metadata: Mapping[str, Any]) -> str:
    metadata = json.dumps(sample_metadata, sort_keys=True).replace("</", "<\\/")
    key = json.dumps(report_key)
    return f"""
const reportKey = {key};
const sampleMetadata = {metadata};
const cards = [...document.querySelectorAll('.route-card')];
const search = document.getElementById('search-filter');
const split = document.getElementById('split-filter');
const steps = document.getElementById('steps-filter');
const status = document.getElementById('status-filter');
const visible = document.getElementById('visible-count');
function storageRead() {{
  try {{ return JSON.parse(localStorage.getItem(reportKey) || '{{}}'); }} catch (_) {{ return {{}}; }}
}}
function storageWrite(value) {{ try {{ localStorage.setItem(reportKey, JSON.stringify(value)); }} catch (_) {{}} }}
let reviews = storageRead();
function applyReviews() {{
  cards.forEach(card => {{
    const routeId = card.dataset.routeId;
    const value = reviews[routeId] || {{status:'unreviewed', note:''}};
    card.dataset.status = value.status || 'unreviewed';
    card.querySelector('.review-status').value = value.status || 'unreviewed';
    card.querySelector('.review-note').value = value.note || '';
  }});
}}
function filterCards() {{
  const query = search.value.trim().toLowerCase(); let count = 0;
  cards.forEach(card => {{
    const show = (!query || card.dataset.search.includes(query)) &&
      (split.value === 'all' || card.dataset.split === split.value) &&
      (steps.value === 'all' || card.dataset.steps === steps.value) &&
      (status.value === 'all' || card.dataset.status === status.value);
    card.classList.toggle('hidden', !show); if (show) count += 1;
  }});
  visible.textContent = `${{count}} / ${{cards.length}} routes`;
}}
document.querySelectorAll('.review-status,.review-note').forEach(control => {{
  control.addEventListener('change', event => {{
    const routeId = event.target.dataset.route;
    const card = cards.find(item => item.dataset.routeId === routeId);
    reviews[routeId] = {{
      status: card.querySelector('.review-status').value,
      note: card.querySelector('.review-note').value
    }};
    card.dataset.status = reviews[routeId].status; storageWrite(reviews); filterCards();
  }});
}});
[search, split, steps, status].forEach(control => control.addEventListener('input', filterCards));
document.getElementById('expand-all').addEventListener('click', () => cards.forEach(card => card.open = true));
document.getElementById('collapse-all').addEventListener('click', () => cards.forEach(card => card.open = false));
document.getElementById('export-review').addEventListener('click', () => {{
  const payload = {{sample:sampleMetadata, reviews:cards.map(card => ({{
    route_id:card.dataset.routeId,
    status:(reviews[card.dataset.routeId] || {{status:'unreviewed'}}).status,
    note:(reviews[card.dataset.routeId] || {{note:''}}).note
  }}))}};
  const blob = new Blob([JSON.stringify(payload, null, 2)], {{type:'application/json'}});
  const link = document.createElement('a'); link.href = URL.createObjectURL(blob);
  link.download = 'route-review.json'; link.click(); URL.revokeObjectURL(link.href);
}});
applyReviews(); filterCards();
"""


def render_route_review_html(
    records: Iterable[Mapping[str, Any]],
    *,
    sample_size: int = 50,
    seed: int = DEFAULT_ROUTE_REVIEW_SEED,
    title: str = "Random route review",
) -> str:
    """Render a random route sample as a self-contained review document."""

    sampled = sample_route_records(records, sample_size=sample_size, seed=seed)
    split_counts = Counter(str(route.get("split") or "unknown") for route in sampled)
    step_counts = Counter(len(route.get("steps") or ()) for route in sampled)
    total_steps = sum(step_counts[key] * key for key in step_counts)
    patent_count = len({str(row.get("patent_id")) for row in sampled})
    abstracted_steps = sum(
        bool(step.get("abstracted_reaction_smiles"))
        for route in sampled
        for step in route.get("steps") or ()
    )
    route_cards = "".join(
        _route_card(route, index)
        for index, route in enumerate(sampled, start=1)
    )
    split_options = "".join(
        f'<option value="{html.escape(value)}">{html.escape(value.title())}</option>'
        for value in sorted(split_counts)
    )
    step_options = "".join(
        f'<option value="{value}">{value} steps</option>'
        for value in sorted(step_counts)
    )
    metadata = {
        "seed": seed,
        "sample_size": sample_size,
        "route_ids": [str(route.get("route_id") or "") for route in sampled],
    }
    report_key = f"route-review:{seed}:{sample_size}"
    return (
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body>"
        f"<header><h1>{html.escape(title)}</h1>"
        f"<p>Simple-random sample of {sample_size} curated routes · seed {seed}</p></header>"
        '<main><section class="summary" aria-label="Sample summary">'
        f'<div><span>Routes</span><strong>{sample_size}</strong></div>'
        f'<div><span>Recorded steps</span><strong>{total_steps}</strong></div>'
        f'<div><span>Abstracted steps</span><strong>{abstracted_steps}</strong></div>'
        f'<div><span>Patents</span><strong>{patent_count}</strong></div>'
        '</section><section class="toolbar" aria-label="Review controls">'
        '<label>Search<input id="search-filter" type="search" placeholder="Route, patent or target"></label>'
        f'<label>Split<select id="split-filter"><option value="all">All splits</option>{split_options}</select></label>'
        f'<label>Length<select id="steps-filter"><option value="all">All lengths</option>{step_options}</select></label>'
        '<label>Status<select id="status-filter"><option value="all">All statuses</option>'
        '<option value="unreviewed">Unreviewed</option><option value="accept">Accept</option>'
        '<option value="question">Questionable</option><option value="reject">Reject</option></select></label>'
        '<button id="expand-all" type="button">Expand all</button>'
        '<button id="collapse-all" type="button">Collapse all</button>'
        '<button id="export-review" type="button">Export review JSON</button>'
        '<span id="visible-count"></span></section>'
        f'<section id="routes">{route_cards}</section></main>'
        f'<script>{_script(report_key, metadata)}</script></body></html>'
    )


def write_route_review_html(
    source_routes: str | Path,
    output_html: str | Path,
    *,
    sample_size: int = 50,
    seed: int = DEFAULT_ROUTE_REVIEW_SEED,
    title: str = "Random route review",
) -> dict[str, Any]:
    """Load, sample, and write a self-contained visual route review."""

    records = load_curated_routes(source_routes)
    document = render_route_review_html(
        records,
        sample_size=sample_size,
        seed=seed,
        title=title,
    )
    output = Path(output_html)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8", newline="\n")
    sampled = sample_route_records(records, sample_size=sample_size, seed=seed)
    return {
        "source_routes": str(Path(source_routes).resolve()),
        "output_html": str(output.resolve()),
        "sample_size": sample_size,
        "seed": seed,
        "route_ids": [str(route.get("route_id") or "") for route in sampled],
        "html_bytes": output.stat().st_size,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render a random route review")
    parser.add_argument("source_routes")
    parser.add_argument("output_html")
    parser.add_argument("--sample-size", type=int, default=50)
    parser.add_argument("--seed", type=int, default=DEFAULT_ROUTE_REVIEW_SEED)
    parser.add_argument("--title", default="Random route review")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the route review renderer."""

    arguments = _parser().parse_args(argv)
    report = write_route_review_html(
        arguments.source_routes,
        arguments.output_html,
        sample_size=arguments.sample_size,
        seed=arguments.seed,
        title=arguments.title,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DEFAULT_ROUTE_REVIEW_SEED",
    "load_curated_routes",
    "render_route_review_html",
    "sample_route_records",
    "write_route_review_html",
]
