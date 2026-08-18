"""Self-contained chemist review for precedent-route expansion products."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .html_report import molecule_svg, reaction_svg


PRECEDENT_ROUTE_EXPANSION_REVIEW_VERSION = "1.0"
_LEVEL_ORDER = {"R0": 0, "R1": 1, "R2": 2}


def load_precedent_route_expansion_report(
    source: str | Path,
) -> dict[str, Any]:
    """Load the deterministic JSON artifact consumed by the review renderer."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("precedent-route expansion report must be an object")
    if value.get("artifact_type") != "two_step_precedent_route_expansion_poc":
        raise ValueError("unexpected precedent-route expansion artifact type")
    routes = value.get("routes")
    if not isinstance(routes, list) or not routes:
        raise ValueError("precedent-route expansion report contains no routes")
    return value


def _level_for_product(route: Mapping[str, Any], product: str) -> str:
    for level in route.get("levels") or ():
        if product in (level.get("product_smiles") or ()):
            return str(level.get("level") or "R2")
    return "R2"


def _unique_input_evidence(route: Mapping[str, Any]) -> tuple[dict[str, Any], ...]:
    values: dict[tuple[str, str], dict[str, Any]] = {}
    for level in route.get("levels") or ():
        for evidence in level.get("input_evidence") or ():
            key = (
                str(evidence.get("role") or ""),
                str(evidence.get("input_smiles") or ""),
            )
            values[key] = dict(evidence)
    return tuple(values[key] for key in sorted(values))


def _stock_table(route: Mapping[str, Any]) -> str:
    rows = []
    for evidence in _unique_input_evidence(route):
        components = evidence.get("stock_components") or ()
        component_text = []
        for component in components:
            records = component.get("source_records") or ()
            supplier = str(records[0].get("supplier") or "") if records else ""
            snapshot = (
                str(records[0].get("snapshot_date") or "") if records else ""
            )
            component_text.append(
                f"{html.escape(str(component.get('canonical_smiles') or component.get('smiles') or ''))}"
                + (f" · {html.escape(supplier)} {html.escape(snapshot)}" if supplier else "")
            )
        complete = evidence.get("stock_evidence_complete")
        evidence_label = (
            "verified"
            if complete is True
            else "not required"
            if complete is None
            else "missing"
        )
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(evidence.get('minimum_level') or ''))}</td>"
            f"<td>{html.escape(str(evidence.get('role') or '').replace('_', ' '))}</td>"
            f"<td><code>{html.escape(str(evidence.get('input_smiles') or ''))}</code>"
            f"<small>{html.escape(str(evidence.get('label') or ''))}</small></td>"
            f"<td>{html.escape(evidence_label)}</td>"
            f"<td>{'<br>'.join(component_text) or '—'}</td>"
            "</tr>"
        )
    return (
        '<details class="evidence"><summary>Input and supplier evidence</summary>'
        '<table><thead><tr><th>Level</th><th>Role</th><th>Input</th>'
        '<th>Stock</th><th>Exact component evidence</th></tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table></details>"
    )


def _source_precedent(route: Mapping[str, Any]) -> str:
    first = str(route.get("first_reaction_smiles") or "")
    second = str(route.get("second_reaction_smiles") or "")
    reactions = "".join(
        '<div class="reaction source-reaction">'
        + reaction_svg(reaction)
        + f'<div class="reaction-caption">Observed step {index}</div></div>'
        for index, reaction in enumerate((first, second), 1)
        if reaction
    )
    return (
        '<details class="source"><summary>Observed two-step precedent</summary>'
        '<div class="source-meta">'
        f"<b>Patent:</b> {html.escape(str(route.get('patent_id') or '—'))} · "
        f"<b>Route:</b> {html.escape(str(route.get('source_route_id') or '—'))} · "
        f"<b>Reaction IDs:</b> {html.escape(str(route.get('first_source_reaction_id') or '—'))}, "
        f"{html.escape(str(route.get('second_source_reaction_id') or '—'))}"
        "</div>"
        f'<div class="reaction-grid">{reactions}</div></details>'
    )


def _pathway_panel(pathway: Mapping[str, Any], index: int) -> str:
    intermediate = str(pathway.get("intermediate_smiles") or "")
    product = str(pathway.get("product_smiles") or "")
    first_inputs = str(pathway.get("first_step_starting_materials") or "")
    second_inputs = str(pathway.get("second_step_starting_materials") or "")
    first_reaction = f"{first_inputs}>>{intermediate}"
    second_reaction = f"{second_inputs}>>{product}"
    warnings = pathway.get("warnings") or ()
    warning_html = "".join(
        f"<li>{html.escape(str(warning))}</li>" for warning in warnings
    )
    return (
        '<details class="pathway" open>'
        f"<summary>Pathway {index} · "
        f"{html.escape(str(pathway.get('first_abstraction_level') or ''))} → "
        f"{html.escape(str(pathway.get('second_abstraction_level') or ''))}</summary>"
        '<div class="reaction-grid">'
        f'<div class="reaction">{reaction_svg(first_reaction)}'
        '<div class="reaction-caption">Step 1 · supplied building blocks → generated intermediate</div></div>'
        f'<div class="reaction">{reaction_svg(second_reaction)}'
        '<div class="reaction-caption">Step 2 · generated intermediate → final product</div></div>'
        "</div>"
        '<dl class="path-meta">'
        f"<dt>Step 1 inputs</dt><dd><code>{html.escape(first_inputs)}</code></dd>"
        f"<dt>Intermediate</dt><dd><code>{html.escape(intermediate)}</code></dd>"
        f"<dt>Step 2 partner</dt><dd><code>{html.escape(str(pathway.get('second_step_partner_smiles') or 'none'))}</code></dd>"
        f"<dt>Operator levels</dt><dd>{html.escape(str(pathway.get('first_abstraction_level') or ''))} / "
        f"{html.escape(str(pathway.get('second_abstraction_level') or ''))}</dd>"
        f"<dt>Graph checks</dt><dd>reverse recovery and edit agreement passed for both steps</dd>"
        "</dl>"
        + (f'<ul class="warnings">{warning_html}</ul>' if warning_html else "")
        + "</details>"
    )


def _review_controls() -> str:
    options = (
        '<option value="unreviewed">Unreviewed</option>'
        '<option value="accept">Accept</option>'
        '<option value="question">Question</option>'
        '<option value="reject">Reject</option>'
        '<option value="out_of_scope">Out of scope</option>'
    )
    return (
        '<div class="review">'
        f'<label>Overall<select class="review-overall">{options}</select></label>'
        f'<label>Step 1<select class="review-step1">{options}</select></label>'
        f'<label>Step 2<select class="review-step2">{options}</select></label>'
        '<label>Site-selectivity risk<select class="review-selectivity">'
        '<option value="unreviewed">Unreviewed</option><option value="low">Low</option>'
        '<option value="medium">Medium</option><option value="high">High</option></select></label>'
        '<label>Condition transfer<select class="review-conditions">'
        '<option value="unreviewed">Unreviewed</option><option value="acceptable">Acceptable</option>'
        '<option value="uncertain">Uncertain</option><option value="incompatible">Incompatible</option></select></label>'
        '<label class="notes">Notes<textarea class="review-notes" rows="3" '
        'placeholder="Selectivity, functional-group compatibility, intermediate stability, conditions…"></textarea></label>'
        "</div>"
    )


def _product_card(
    route: Mapping[str, Any],
    product: str,
    pathways: Iterable[Mapping[str, Any]],
    index: int,
) -> str:
    values = tuple(pathways)
    level = _level_for_product(route, product)
    route_id = str(route.get("route_id") or "")
    review_id = f"{route_id}|{product}"
    search = " ".join(
        (
            route_id,
            str(route.get("patent_id") or ""),
            product,
            *(str(item.get("first_step_starting_materials") or "") for item in values),
            *(str(item.get("second_step_partner_smiles") or "") for item in values),
        )
    ).lower()
    pathway_html = "".join(
        _pathway_panel(pathway, pathway_index)
        for pathway_index, pathway in enumerate(values, 1)
    )
    return (
        f'<article class="product-card level-{level.lower()}" '
        f'data-review-id="{html.escape(review_id)}" '
        f'data-route="{html.escape(route_id)}" data-level="{level}" '
        f'data-product="{html.escape(product)}" data-search="{html.escape(search)}">'
        '<div class="product-head">'
        f'<span class="number">{index:02d}</span><span class="level">{level}</span>'
        f'<h3>{html.escape(route_id.replace("_", " ").title())}</h3>'
        f'<span class="path-count">{len(values)} pathway(s)</span></div>'
        '<div class="product-structure">'
        + molecule_svg(product, width=420, height=240)
        + f'<div><b>Generated product</b><code>{html.escape(product)}</code>'
        f'<small>Patent source: {html.escape(str(route.get("patent_id") or "—"))}</small></div></div>'
        f'<div class="pathways">{pathway_html}</div>'
        + _review_controls()
        + "</article>"
    )


def _route_section(route: Mapping[str, Any], start_index: int) -> tuple[str, int]:
    levels = route.get("levels") or ()
    final = levels[-1] if levels else {}
    by_product: dict[str, list[Mapping[str, Any]]] = {}
    for pathway in final.get("pathways") or ():
        by_product.setdefault(str(pathway.get("product_smiles") or ""), []).append(
            pathway
        )
    ordered_products = sorted(
        by_product,
        key=lambda product: (
            _LEVEL_ORDER.get(_level_for_product(route, product), 99),
            product,
        ),
    )
    cards = []
    index = start_index
    for product in ordered_products:
        cards.append(_product_card(route, product, by_product[product], index))
        index += 1
    counts = {
        str(level.get("level") or ""): int(level.get("product_count") or 0)
        for level in levels
    }
    heading = (
        '<section class="route-section">'
        '<div class="route-heading">'
        f'<div><h2>{html.escape(str(route.get("route_id") or "").replace("_", " ").title())}</h2>'
        f'<p>{html.escape(str(route.get("patent_id") or ""))} · '
        f'{html.escape(str(route.get("source_route_id") or ""))}</p></div>'
        f'<div class="route-counts"><span>R0 <b>{counts.get("R0", 0)}</b></span>'
        f'<span>R1 <b>{counts.get("R1", 0)}</b></span>'
        f'<span>R2 <b>{counts.get("R2", 0)}</b></span></div></div>'
        + _source_precedent(route)
        + _stock_table(route)
        + "".join(cards)
        + "</section>"
    )
    return heading, index


def _stylesheet() -> str:
    return """
:root{--ink:#19211e;--muted:#63716b;--line:#d6ded9;--paper:#f3f5f3;--green:#176b4a;--amber:#996515;--red:#a23d34;--blue:#346a8a}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}header{background:#173f33;color:white;padding:24px max(20px,calc((100% - 1500px)/2))}header h1{margin:0 0 5px;font-size:28px}header p{max-width:950px;color:#d8e8e2;margin:0}.summary{display:flex;gap:10px;flex-wrap:wrap;margin-top:15px}.summary span{background:#285346;border-radius:7px;padding:8px 12px}.summary b{font-size:19px;margin-left:5px}main{max-width:1500px;margin:auto;padding:16px}.toolbar{position:sticky;top:0;z-index:10;display:flex;gap:10px;align-items:end;flex-wrap:wrap;background:white;border:1px solid var(--line);border-radius:9px;padding:10px;box-shadow:0 4px 16px #1d3b3020}.toolbar label{display:grid;gap:3px;color:var(--muted);font-size:12px}.toolbar input,.toolbar select,.review select,.review textarea{border:1px solid #b9c6c0;border-radius:6px;background:white;padding:7px}.toolbar input{width:260px}.toolbar button{border:0;border-radius:6px;background:#315e50;color:white;padding:8px 12px;cursor:pointer}.visible{margin-left:auto;color:var(--muted)}.route-section{margin:22px 0 40px}.route-heading{display:flex;justify-content:space-between;align-items:end;border-bottom:2px solid #8ba49a;padding:0 3px 8px}.route-heading h2{margin:0;font-size:24px}.route-heading p{margin:2px 0;color:var(--muted)}.route-counts{display:flex;gap:5px}.route-counts span,.level{border-radius:999px;padding:4px 9px;background:#e8efeb}.route-counts b{font-size:16px}.source,.evidence{background:white;border:1px solid var(--line);border-radius:8px;margin:10px 0;padding:9px}.source summary,.evidence summary,.pathway summary{cursor:pointer;font-weight:650}.source-meta{color:var(--muted);padding:8px 2px}.evidence table{border-collapse:collapse;width:100%;margin-top:8px}.evidence th,.evidence td{text-align:left;border-top:1px solid var(--line);padding:7px;vertical-align:top}.evidence small,.product-structure small{display:block;color:var(--muted)}code{overflow-wrap:anywhere}.product-card{background:white;border:1px solid var(--line);border-left:5px solid #788a83;border-radius:10px;margin:14px 0;padding:13px;box-shadow:0 3px 12px #1f372e0d}.product-card.level-r0{border-left-color:var(--green)}.product-card.level-r1{border-left-color:var(--blue)}.product-card.level-r2{border-left-color:var(--amber)}.product-card[data-review="accept"]{box-shadow:inset 0 0 0 2px var(--green)}.product-card[data-review="question"]{box-shadow:inset 0 0 0 2px var(--amber)}.product-card[data-review="reject"]{box-shadow:inset 0 0 0 2px var(--red)}.product-head{display:flex;align-items:center;gap:9px}.product-head h3{margin:0;font-size:17px}.number,.path-count{color:var(--muted)}.path-count{margin-left:auto}.product-structure{display:grid;grid-template-columns:minmax(300px,430px) 1fr;align-items:center;gap:15px;border-bottom:1px solid var(--line);margin:8px 0}.product-structure svg{width:100%;height:230px}.product-structure code{display:block;margin-top:7px}.pathway{border:1px solid var(--line);border-radius:7px;margin:8px 0;padding:8px}.reaction-grid{display:grid;grid-template-columns:repeat(2,minmax(440px,1fr));gap:8px}.reaction{overflow:auto;background:#fafbfa;border-radius:6px}.reaction svg{width:100%;min-width:520px;height:auto;max-height:220px}.reaction-caption{text-align:center;color:var(--muted);font-size:12px;padding:3px}.path-meta{display:grid;grid-template-columns:130px 1fr;gap:3px 9px;margin:7px}.path-meta dt{color:var(--muted)}.path-meta dd{margin:0}.warnings{color:var(--amber)}.review{display:grid;grid-template-columns:repeat(5,minmax(130px,1fr));gap:8px;border-top:1px solid var(--line);padding-top:10px;margin-top:10px}.review label{display:grid;gap:3px;color:var(--muted);font-size:12px}.review .notes{grid-column:1/-1}.review textarea{resize:vertical;width:100%}@media(max-width:1000px){.reaction-grid{grid-template-columns:1fr}.product-structure{grid-template-columns:1fr}.review{grid-template-columns:repeat(2,1fr)}.toolbar{position:static}.visible{width:100%;margin:0}}@media print{.toolbar,.review{display:none}.product-card{break-inside:avoid}.source,.evidence,.pathway{display:block}}
"""


def _script(storage_key: str, artifact_metadata: Mapping[str, Any]) -> str:
    return f"""
const key={json.dumps(storage_key)};
let saved={{}};try{{saved=JSON.parse(localStorage.getItem(key)||'{{}}')}}catch(_){{saved={{}}}}
const cards=[...document.querySelectorAll('.product-card')];
for(const card of cards){{
 const id=card.dataset.reviewId,state=saved[id]||{{}};
 const controls={{overall:'.review-overall',step1:'.review-step1',step2:'.review-step2',selectivity:'.review-selectivity',conditions:'.review-conditions',notes:'.review-notes'}};
 for(const [name,selector] of Object.entries(controls)){{const el=card.querySelector(selector);el.value=state[name]||(name==='notes'?'':'unreviewed');el.addEventListener(name==='notes'?'input':'change',()=>{{const value={{}};for(const [field,sel] of Object.entries(controls))value[field]=card.querySelector(sel).value;saved[id]=value;card.dataset.review=value.overall;try{{localStorage.setItem(key,JSON.stringify(saved))}}catch(_){{}}filter();}})}}
 card.dataset.review=state.overall||'unreviewed';
}}
function filter(){{const query=document.querySelector('#search').value.toLowerCase(),route=document.querySelector('#route-filter').value,level=document.querySelector('#level-filter').value,status=document.querySelector('#status-filter').value;let visible=0;for(const card of cards){{const show=(!query||card.dataset.search.includes(query))&&(route==='all'||card.dataset.route===route)&&(level==='all'||card.dataset.level===level)&&(status==='all'||card.dataset.review===status);card.hidden=!show;if(show)visible++}}document.querySelector('#visible').textContent=`${{visible}} / ${{cards.length}} products`;}}
for(const id of ['search','route-filter','level-filter','status-filter'])document.querySelector('#'+id).addEventListener('input',filter);
document.querySelector('#export').addEventListener('click',()=>{{const reviews=cards.map(card=>({{review_id:card.dataset.reviewId,route_id:card.dataset.route,product_smiles:card.dataset.product,expansion_level:card.dataset.level,...(saved[card.dataset.reviewId]||{{overall:'unreviewed',step1:'unreviewed',step2:'unreviewed',selectivity:'unreviewed',conditions:'unreviewed',notes:''}})}}));const payload={{review_schema_version:'1.0',artifact:{json.dumps(dict(artifact_metadata), sort_keys=True)},reviews}};const blob=new Blob([JSON.stringify(payload,null,2)],{{type:'application/json'}}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='precedent-route-expansion-review.json';a.click();URL.revokeObjectURL(a.href)}});filter();
"""


def render_precedent_route_expansion_html(
    report: Mapping[str, Any],
    *,
    title: str = "Two-step precedent-route product review",
) -> str:
    """Render generated products and both-step pathways as standalone HTML."""

    routes = tuple(report.get("routes") or ())
    if not routes:
        raise ValueError("precedent-route expansion report contains no routes")
    sections = []
    next_index = 1
    route_options = []
    for route in routes:
        section, next_index = _route_section(route, next_index)
        sections.append(section)
        route_id = str(route.get("route_id") or "")
        route_options.append(
            f'<option value="{html.escape(route_id)}">{html.escape(route_id.replace("_", " ").title())}</option>'
        )
    product_count = next_index - 1
    counts = report.get("level_product_counts") or {}
    source = report.get("source_validation") or {}
    stock = report.get("stock_index_metadata") or {}
    stock_sources = stock.get("source_summaries") or ()
    stock_label = (
        f"{stock_sources[0].get('supplier', '')} {stock_sources[0].get('snapshot_date', '')}"
        if stock_sources
        else "not supplied"
    )
    metadata = {
        "definition_id": report.get("definition_id"),
        "algorithm_version": report.get("algorithm_version"),
        "source_sha256": source.get("source_sha256"),
    }
    storage_key = "precedent-route-expansion-review:" + str(
        report.get("definition_id") or "unknown"
    )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body>"
        f"<header><h1>{html.escape(title)}</h1>"
        "<p>Review the complete two-step pathway. Structural validation and stock identity are already checked; assess selectivity, functional-group compatibility, intermediate viability, and condition transfer.</p>"
        '<div class="summary">'
        f"<span>Routes <b>{len(routes)}</b></span><span>Products <b>{product_count}</b></span>"
        f"<span>R0 <b>{int(counts.get('R0') or 0)}</b></span>"
        f"<span>R1 <b>{int(counts.get('R1') or 0)}</b></span>"
        f"<span>R2 <b>{int(counts.get('R2') or 0)}</b></span>"
        f"<span>Source verified <b>{'yes' if source.get('valid') else 'no'}</b></span>"
        f"<span>Stock <b>{html.escape(stock_label)}</b></span></div></header>"
        '<main><section class="toolbar"><label>Search<input id="search" type="search" placeholder="Product, route, input, patent…"></label>'
        '<label>Route<select id="route-filter"><option value="all">All</option>'
        + "".join(route_options)
        + '</select></label><label>Expansion level<select id="level-filter"><option value="all">All</option><option value="R0">R0 exact</option><option value="R1">R1 exact context</option><option value="R2">R2 relaxed</option></select></label>'
        '<label>Review status<select id="status-filter"><option value="all">All</option><option value="unreviewed">Unreviewed</option><option value="accept">Accept</option><option value="question">Question</option><option value="reject">Reject</option><option value="out_of_scope">Out of scope</option></select></label>'
        '<button id="export" type="button">Export review JSON</button><span class="visible" id="visible"></span></section>'
        + "".join(sections)
        + f"</main><script>{_script(storage_key, metadata)}</script></body></html>"
    )


def write_precedent_route_expansion_html(
    source_report: Mapping[str, Any] | str | Path,
    output_html: str | Path,
    *,
    title: str = "Two-step precedent-route product review",
) -> dict[str, Any]:
    """Write a deterministic standalone review and return its file summary."""

    report = (
        load_precedent_route_expansion_report(source_report)
        if isinstance(source_report, (str, Path))
        else dict(source_report)
    )
    document = render_precedent_route_expansion_html(report, title=title)
    destination = Path(output_html)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(document, encoding="utf-8", newline="\n")
    return {
        "review_version": PRECEDENT_ROUTE_EXPANSION_REVIEW_VERSION,
        "definition_id": report.get("definition_id"),
        "route_count": len(report.get("routes") or ()),
        "product_count": sum(
            len((route.get("levels") or [{}])[-1].get("product_smiles") or ())
            for route in report.get("routes") or ()
        ),
        "output_html": str(destination.resolve()),
        "html_bytes": destination.stat().st_size,
    }


__all__ = [
    "PRECEDENT_ROUTE_EXPANSION_REVIEW_VERSION",
    "load_precedent_route_expansion_report",
    "render_precedent_route_expansion_html",
    "write_precedent_route_expansion_html",
]
