"""Self-contained HTML review for retrosynthesis comparison artifacts."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from visualization import (
    render_molecule_image_bytes,
    render_reaction_image_bytes,
)


DEFAULT_METHODS = (
    "baseline",
    "core_l1_context",
    "ensemble_baseline_l1_context",
)

METHOD_LABELS = {
    "baseline": "RDChiral baseline",
    "core_l1_context": "Core L1 + context",
    "core_l1_neutral": "Core L1 neutral",
    "core_l2_context": "Core L2 + context",
    "core_l2_neutral": "Core L2 neutral",
    "core_context": "Core L1/L2 + context",
    "core_site_diverse": "Core + site-diverse ranking",
    "core_neutral": "Core L1/L2 neutral",
    "ensemble_baseline_l1_context": "Baseline-first L1 ensemble",
    "ensemble_baseline_core_context": "Baseline-first generic-core ensemble",
    "ensemble_baseline_core_site_diverse": (
        "Baseline-first site-diverse ensemble"
    ),
    "supported_l1_l2": "Seven-archetype L1/L2",
    "data_l2": "Data-derived L2",
    "data_l2_l1": "Data-derived L2/L1",
    "data_ladder": "Data-derived L2/L1/L0",
}


def _inline_svg(value: str) -> str:
    start = value.find("<svg")
    return value[start:] if start >= 0 else value


def _placeholder(message: str) -> str:
    return f'<div class="render-error">{html.escape(message)}</div>'


def molecule_svg(smiles: str, *, width: int = 320, height: int = 220) -> str:
    """Render one molecule as an inline SVG or a visible error placeholder."""

    try:
        drawing = render_molecule_image_bytes(
            smiles,
            size=(width, height),
            image_format="svg",
        )
        return _inline_svg(drawing.decode("utf-8"))
    except (RuntimeError, ValueError, UnicodeDecodeError):
        return _placeholder("Molecule drawing failed")


def reaction_svg(
    reaction_smiles: str,
    *,
    sub_image_size: tuple[int, int] = (260, 180),
) -> str:
    """Render one reaction SMILES as an inline SVG or an error placeholder."""

    try:
        drawing = render_reaction_image_bytes(
            reaction_smiles,
            size=(sub_image_size[0] * 3, sub_image_size[1]),
            image_format="svg",
        )
        return _inline_svg(drawing.decode("utf-8"))
    except (RuntimeError, ValueError, UnicodeDecodeError):
        return _placeholder("Reaction drawing failed")


def _metric(value: Any) -> str:
    try:
        return f"{float(value):.3f}"
    except (TypeError, ValueError):
        return "—"


def _method_label(method: str) -> str:
    return METHOD_LABELS.get(method, method.replace("_", " ").title())


def _metrics_table(data: dict[str, Any], methods: Sequence[str]) -> str:
    metrics = data.get("metrics") or {}
    has_site_metrics = any(
        "top1_site_recall" in (metrics.get(method) or {})
        for method in methods
    )
    has_operator_metrics = any(
        "top1_operator_recall" in (metrics.get(method) or {})
        for method in methods
    )
    rows = []
    for method in methods:
        values = metrics.get(method) or {}
        if has_operator_metrics:
            recall_cells = (
                f"<td>{_metric(values.get('top1_exact_precursor_recall'))}</td>"
                f"<td>{_metric(values.get('top25_exact_precursor_recall'))}</td>"
                f"<td>{_metric(values.get('top1_synthon_recall'))}</td>"
                f"<td>{_metric(values.get('top25_synthon_recall'))}</td>"
                f"<td>{_metric(values.get('top1_operator_recall'))}</td>"
                f"<td>{_metric(values.get('top25_operator_recall'))}</td>"
                f"<td>{_metric(values.get('top1_site_recall'))}</td>"
                f"<td>{_metric(values.get('top25_site_recall'))}</td>"
            )
        elif has_site_metrics:
            recall_cells = (
                f"<td>{_metric(values.get('top1_exact_precursor_recall'))}</td>"
                f"<td>{_metric(values.get('top10_exact_precursor_recall'))}</td>"
                f"<td>{_metric(values.get('top1_site_recall'))}</td>"
                f"<td>{_metric(values.get('top5_site_recall'))}</td>"
                f"<td>{_metric(values.get('top10_site_recall'))}</td>"
            )
        else:
            recall_cells = (
                f"<td>{_metric(values.get('top1_exact_precursor_recall'))}</td>"
                f"<td>{_metric(values.get('top5_exact_precursor_recall'))}</td>"
                f"<td>{_metric(values.get('top10_exact_precursor_recall'))}</td>"
            )
        rows.append(
            "<tr>"
            f"<th>{html.escape(_method_label(method))}</th>"
            f"<td>{_metric(values.get('target_coverage'))}</td>"
            f"{recall_cells}"
            f"<td>{_metric(values.get('valid_candidate_fraction'))}</td>"
            f"<td>{_metric(values.get('mean_candidates_per_target'))}</td>"
            "</tr>"
        )
    recall_headers = (
        "<th>Exact @1</th><th>Exact @25</th><th>Synthon @1</th>"
        "<th>Synthon @25</th><th>Operator @1</th><th>Operator @25</th>"
        "<th>Site @1</th><th>Site @25</th>"
        if has_operator_metrics
        else
        "<th>Exact @1</th><th>Exact @10</th><th>Site @1</th>"
        "<th>Site @5</th><th>Site @10</th>"
        if has_site_metrics
        else "<th>Top 1</th><th>Top 5</th><th>Top 10</th>"
    )
    return (
        '<table class="metrics"><thead><tr><th>Method</th><th>Coverage</th>'
        f"{recall_headers}<th>Valid</th>"
        f"<th>Candidates/target</th></tr></thead><tbody>{''.join(rows)}</tbody></table>"
    )


def _rank_badge(label: str, rank: Any, css_class: str) -> str:
    value = "miss" if rank is None else f"rank {rank}"
    return (
        f'<span class="badge {css_class}">{html.escape(label)}: '
        f"{html.escape(value)}</span>"
    )


def _candidate_card(
    precursor_smiles: str,
    target_smiles: str,
    rank: int,
    exact_rank: int | None,
    site_rank: int | None,
    operator_rank: int | None = None,
    synthon_rank: int | None = None,
) -> str:
    reaction_smiles = f"{precursor_smiles}>>{target_smiles}"
    badges = [f'<span class="badge rank">#{rank}</span>']
    if exact_rank == rank:
        badges.append('<span class="badge exact">exact precursor</span>')
    if site_rank == rank:
        badges.append('<span class="badge center">correct product site</span>')
    if operator_rank == rank:
        badges.append('<span class="badge operator">correct operator</span>')
    if synthon_rank == rank:
        badges.append('<span class="badge synthon">equivalent synthons</span>')
    return (
        '<article class="candidate">'
        f'<div class="candidate-badges">{"".join(badges)}</div>'
        f'<div class="reaction-svg">{reaction_svg(reaction_smiles)}</div>'
        f'<code title="{html.escape(reaction_smiles)}">'
        f"{html.escape(reaction_smiles)}</code>"
        "</article>"
    )


def _method_section(
    target: dict[str, Any],
    method: str,
    top_k: int,
) -> str:
    values = (target.get("methods") or {}).get(method) or {}
    precursors = tuple(values.get("precursor_smiles") or ())[:top_k]
    exact_rank = values.get("exact_precursor_rank")
    site_rank = values.get("site_rank", values.get("center_rank"))
    operator_rank = values.get("operator_rank")
    synthon_rank = values.get("synthon_rank")
    candidates = "".join(
        _candidate_card(
            str(precursor),
            str(target.get("target_smiles") or ""),
            rank,
            exact_rank,
            site_rank,
            operator_rank,
            synthon_rank,
        )
        for rank, precursor in enumerate(precursors, start=1)
    )
    if not candidates:
        candidates = '<p class="empty">No candidate generated.</p>'
    return (
        '<section class="method-section">'
        '<div class="method-title">'
        f"<h4>{html.escape(_method_label(method))}</h4>"
        f'{_rank_badge("exact", exact_rank, "exact-rank")}'
        f'{_rank_badge("site", site_rank, "center-rank")}'
        f'{_rank_badge("operator", operator_rank, "operator-rank")}'
        f'{_rank_badge("synthon", synthon_rank, "synthon-rank")}'
        "</div>"
        f'<div class="candidate-grid">{candidates}</div>'
        "</section>"
    )


def _target_card(
    target: dict[str, Any],
    methods: Sequence[str],
    top_k: int,
    index: int,
) -> str:
    target_smiles = str(target.get("target_smiles") or "")
    expected = str(target.get("expected_precursor_smiles") or "")
    bond_kind = str(target.get("bond_kind") or "unknown")
    reaction_id = str(target.get("reaction_id") or f"target-{index}")
    baseline_rank = ((target.get("methods") or {}).get("baseline") or {}).get(
        "exact_precursor_rank"
    )
    target_methods = target.get("methods") or {}
    ensemble_method = next(
        (method for method in methods if method.startswith("ensemble_")),
        methods[-1] if methods else "baseline",
    )
    ensemble_rank = (target_methods.get(ensemble_method) or {}).get(
        "exact_precursor_rank"
    )
    recovery = (
        "rescue"
        if baseline_rank is None and ensemble_rank is not None
        else "recovered"
        if ensemble_rank is not None
        else "missed"
    )
    summary_badges = (
        f'<span class="badge bond">{html.escape(bond_kind)}</span>'
        f'<span class="badge recovery {recovery}">{recovery}</span>'
        f'{_rank_badge("ensemble exact", ensemble_rank, "exact-rank")}'
    )
    expected_reaction = f"{expected}>>{target_smiles}"
    method_sections = "".join(
        _method_section(target, method, top_k) for method in methods
    )
    search_value = " ".join((reaction_id, target_smiles, expected)).lower()
    return (
        f'<details class="target-card" data-bond="{html.escape(bond_kind)}" '
        f'data-recovery="{recovery}" data-search="{html.escape(search_value)}"'
        f"{' open' if index == 1 else ''}>"
        "<summary>"
        f'<span class="target-index">{index}</span>'
        f'<strong>{html.escape(reaction_id)}</strong>{summary_badges}'
        f'<code>{html.escape(target_smiles)}</code>'
        "</summary>"
        '<div class="target-body">'
        '<section class="query-panel">'
        '<div><h3>Query compound</h3>'
        f'<div class="query-svg">{molecule_svg(target_smiles)}</div>'
        f'<code>{html.escape(target_smiles)}</code></div>'
        '<div><h3>Recorded reaction</h3>'
        f'<div class="expected-svg">{reaction_svg(expected_reaction)}</div>'
        f'<code>{html.escape(expected_reaction)}</code></div>'
        "</section>"
        f"{method_sections}</div></details>"
    )


def _stylesheet() -> str:
    return """
:root { color-scheme: light; --ink:#17212b; --muted:#617080; --line:#dbe2e8;
  --panel:#fff; --wash:#f4f7f9; --accent:#176b87; --good:#176b4d;
  --warn:#9a6700; --bad:#a02d2d; }
* { box-sizing:border-box; }
body { margin:0; background:var(--wash); color:var(--ink); font:14px/1.45
  Inter,Segoe UI,Arial,sans-serif; }
header { background:linear-gradient(120deg,#12394b,#176b87); color:#fff;
  padding:28px max(24px,calc((100vw - 1500px)/2)); }
header h1 { margin:0 0 7px; font-size:28px; } header p { margin:0; opacity:.86; }
main { max-width:1500px; margin:auto; padding:22px; }
.summary-panel,.toolbar,.target-card { background:var(--panel); border:1px solid var(--line);
  border-radius:12px; box-shadow:0 2px 10px #12202b0c; }
.summary-panel { padding:18px; margin-bottom:16px; overflow:auto; }
.metrics { width:100%; border-collapse:collapse; min-width:760px; }
.metrics th,.metrics td { padding:9px 11px; border-bottom:1px solid var(--line);
  text-align:right; } .metrics th:first-child { text-align:left; }
.toolbar { position:sticky; top:0; z-index:4; display:flex; gap:12px; align-items:end;
  padding:12px 16px; margin-bottom:16px; }
.toolbar label { display:grid; gap:4px; color:var(--muted); font-size:12px; }
.toolbar input,.toolbar select { min-width:170px; padding:8px 10px; border:1px solid
  #bcc8d1; border-radius:7px; background:#fff; color:var(--ink); }
#visible-count { margin-left:auto; color:var(--muted); padding-bottom:7px; }
.target-card { margin-bottom:13px; overflow:hidden; }
.target-card>summary { cursor:pointer; display:flex; flex-wrap:wrap; gap:9px;
  align-items:center; padding:13px 16px; background:#fbfcfd; }
.target-card>summary code { flex-basis:100%; margin-left:36px; }
.target-index { display:grid; place-items:center; width:27px; height:27px; border-radius:50%;
  background:#dcecf2; color:#124c62; font-weight:700; }
.target-body { padding:16px; }.query-panel { display:grid;
  grid-template-columns:minmax(280px,.7fr) minmax(450px,1.3fr); gap:16px; }
.query-panel>div,.method-section { border:1px solid var(--line); border-radius:10px;
  padding:12px; overflow:hidden; }.query-panel h3,.method-title h4 { margin:0 0 8px; }
.query-svg svg,.expected-svg svg,.reaction-svg svg { width:100%; height:auto; max-height:260px; }
.method-section { margin-top:14px; }.method-title { display:flex; gap:8px;
  align-items:center; flex-wrap:wrap; }.method-title h4 { margin-right:auto; font-size:17px; }
.candidate-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(330px,1fr)); gap:10px; }
.candidate { position:relative; border:1px solid #e4e9ed; background:#fcfdfd;
  border-radius:8px; padding:9px; overflow:hidden; }.candidate-badges { min-height:25px; }
code { display:block; white-space:normal; overflow-wrap:anywhere; color:#33414d;
  font:12px/1.35 Consolas,monospace; }.badge { display:inline-block; padding:3px 7px;
  margin:1px 3px 1px 0; border-radius:999px; background:#e8edf1; font-size:11px;
  font-weight:650; }.badge.exact,.recovery.recovered { background:#d9f2e7; color:var(--good); }
.badge.center { background:#dcecf7; color:#155477; }.badge.recovery.rescue {
  background:#fff0c7; color:var(--warn); }.badge.recovery.missed { background:#f9dddd;
  color:var(--bad); }.badge.bond { background:#e6e2fa; color:#4b3884; }
.badge.operator { background:#e2edf9; color:#23537a; }
.badge.synthon { background:#e5f2dd; color:#35651e; }
.render-error,.empty { padding:30px; color:var(--bad); text-align:center; }
.hidden { display:none; }
@media (max-width:800px) { .query-panel { grid-template-columns:1fr; }
  .toolbar { position:static; flex-wrap:wrap; } .candidate-grid { grid-template-columns:1fr; } }
"""


def _script() -> str:
    return """
const bond = document.getElementById('bond-filter');
const recovery = document.getElementById('recovery-filter');
const search = document.getElementById('search-filter');
const cards = [...document.querySelectorAll('.target-card')];
const visible = document.getElementById('visible-count');
function filterCards() {
  const query = search.value.trim().toLowerCase(); let count = 0;
  cards.forEach(card => {
    const show = (bond.value === 'all' || card.dataset.bond === bond.value)
      && (recovery.value === 'all' || card.dataset.recovery === recovery.value)
      && (!query || card.dataset.search.includes(query));
    card.classList.toggle('hidden', !show); if (show) count += 1;
  });
  visible.textContent = `${count} of ${cards.length} targets`;
}
[bond,recovery,search].forEach(control => control.addEventListener('input', filterCards));
filterCards();
"""


def render_comparison_html(
    comparison: dict[str, Any],
    *,
    methods: Iterable[str] = DEFAULT_METHODS,
    top_k: int = 5,
    max_targets: int | None = None,
    title: str = "Retrosynthesis comparison review",
) -> str:
    """Return one self-contained, filterable comparison HTML document."""

    if top_k < 1:
        raise ValueError("top-k must be positive")
    if max_targets is not None and max_targets < 1:
        raise ValueError("max targets must be positive")
    available = comparison.get("metrics") or {}
    selected_methods = tuple(dict.fromkeys(methods))
    missing = [method for method in selected_methods if method not in available]
    if missing:
        raise ValueError(f"comparison does not contain methods: {', '.join(missing)}")
    targets = tuple(comparison.get("target_results") or ())[:max_targets]
    bond_options = "".join(
        f'<option value="{html.escape(value)}">{html.escape(value)}</option>'
        for value in sorted(
            {
                str(target.get("bond_kind") or "unclassified")
                for target in targets
            }
        )
    )
    cards = "".join(
        _target_card(target, selected_methods, top_k, index)
        for index, target in enumerate(targets, start=1)
    )
    split = comparison.get("split") or {}
    subtitle = (
        f"{len(targets)} rendered targets; source rows "
        f"{split.get('source_rows', '—')}; held-out groups "
        f"{split.get('held_out_support_groups', '—')}; top {top_k} shown."
    )
    return (
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body>"
        f"<header><h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p></header>"
        '<main><section class="summary-panel"><h2>Method summary</h2>'
        f"{_metrics_table(comparison, selected_methods)}</section>"
        '<section class="toolbar"><label>Transformation<select id="bond-filter">'
        f'<option value="all">All classes</option>{bond_options}</select></label>'
        '<label>Recovery<select id="recovery-filter">'
        '<option value="all">All outcomes</option><option value="rescue">Core rescues</option>'
        '<option value="recovered">Recovered</option><option value="missed">Missed</option>'
        '</select></label><label>Search<input id="search-filter" '
        'placeholder="Reaction ID or SMILES"></label><span id="visible-count"></span></section>'
        f'<section id="targets">{cards}</section></main><script>{_script()}</script>'
        "</body></html>"
    )


def write_comparison_html(
    comparison_source: str | Path,
    output: str | Path,
    *,
    methods: Iterable[str] = DEFAULT_METHODS,
    top_k: int = 5,
    max_targets: int | None = None,
    title: str = "Retrosynthesis comparison review",
) -> dict[str, Any]:
    """Read comparison JSON and write a self-contained HTML review."""

    source_path = Path(comparison_source)
    data = json.loads(source_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("comparison JSON must contain an object")
    selected_methods = tuple(dict.fromkeys(methods))
    document = render_comparison_html(
        data,
        methods=selected_methods,
        top_k=top_k,
        max_targets=max_targets,
        title=title,
    )
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8", newline="\n")
    target_count = min(
        len(data.get("target_results") or ()),
        max_targets if max_targets is not None else 10**12,
    )
    return {
        "output": str(output_path),
        "target_count": target_count,
        "methods": list(selected_methods),
        "top_k": top_k,
        "self_contained": True,
    }


__all__ = [
    "DEFAULT_METHODS",
    "molecule_svg",
    "reaction_svg",
    "render_comparison_html",
    "write_comparison_html",
]
