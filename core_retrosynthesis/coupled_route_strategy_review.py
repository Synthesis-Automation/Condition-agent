"""Self-contained chemist review for coupled route-strategy mining."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Mapping

from .coupled_route_strategy import COUPLED_ROUTE_STRATEGY_REVIEW_VERSION
from .html_report import molecule_svg, reaction_svg


def load_coupled_route_strategy_report(source: str | Path) -> dict[str, Any]:
    """Load and validate one coupled-strategy POC report."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("coupled route-strategy report must be an object")
    if value.get("artifact_type") != "coupled_route_strategy_poc":
        raise ValueError("unexpected coupled route-strategy artifact type")
    if not isinstance(value.get("sample"), list) or not value["sample"]:
        raise ValueError("coupled route-strategy report has no review sample")
    return value


def _text_list(values: Any) -> str:
    return "".join(f"<li><code>{html.escape(str(item))}</code></li>" for item in values)


def _strategy_summary_table(report: Mapping[str, Any]) -> str:
    rows = []
    for item in (report.get("strategy_summaries") or ())[:12]:
        first = "<br>".join(
            html.escape(str(value)) for value in item.get("first_edit_tokens") or ()
        )
        second = "<br>".join(
            html.escape(str(value)) for value in item.get("second_edit_tokens") or ()
        )
        rows.append(
            "<tr>"
            f'<td>{html.escape(str(item.get("dependency_class") or "").replace("_", " "))}</td>'
            f'<td>{int(item.get("patent_count") or 0)}</td>'
            f'<td>{int(item.get("occurrence_count") or 0)}</td>'
            f"<td><code>{first}</code></td><td><code>{second}</code></td>"
            "</tr>"
        )
    if not rows:
        return ""
    return (
        '<details class="top-strategies"><summary>Most recurrent strict typed strategies</summary>'
        '<table><thead><tr><th>Dependency</th><th>Patents</th><th>Occurrences</th>'
        '<th>Step 1 edits</th><th>Step 2 edits</th></tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table></details>'
    )


def _comparison_table(report: Mapping[str, Any]) -> str:
    rows = []
    for item in report.get("baseline_to_dependency_counts") or ():
        rows.append(
            "<tr>"
            f'<td>{html.escape(str(item.get("relationship_class") or "").replace("_", " "))}</td>'
            '<td class="arrow">→</td>'
            f'<td>{html.escape(str(item.get("dependency_class") or "").replace("_", " "))}</td>'
            f'<td>{int(item.get("count") or 0):,}</td>'
            "</tr>"
        )
    if not rows:
        return ""
    return (
        '<details class="comparison" open><summary>Baseline relationship → refined dependency</summary>'
        '<p>The baseline asks whether the two active sites overlap. The refined class asks how step 2 depends on the handle made or transformed by step 1.</p>'
        '<table><thead><tr><th>Baseline v1</th><th></th><th>Refined v2</th><th>Pairs</th></tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table></details>'
    )


def _evidence_panel(item: Mapping[str, Any]) -> str:
    evidence = item.get("evidence") or {}
    distances = [
        "unresolved" if value is None else str(value)
        for value in evidence.get("minimum_distances") or ()
    ]
    rationale = "".join(
        f"<li>{html.escape(str(value))}</li>"
        for value in evidence.get("rationale") or ()
    )
    warnings = "".join(
        f"<li>{html.escape(str(value))}</li>"
        for value in item.get("warnings") or ()
    )
    return (
        '<div class="evidence-grid">'
        f'<div><b>Step-1 active maps</b><code>{html.escape(str(evidence.get("producer_active_maps") or []))}</code></div>'
        f'<div><b>Persistent active maps</b><code>{html.escape(str(evidence.get("producer_persistent_active_maps") or []))}</code></div>'
        f'<div><b>Step-2 intermediate maps</b><code>{html.escape(str(evidence.get("consumer_intermediate_active_maps") or []))}</code></div>'
        f'<div><b>Overlap by lineage</b><code>{html.escape(str(evidence.get("overlap_counts") or []))}</code></div>'
        f'<div><b>Formed-bond overlap</b><code>{html.escape(str(evidence.get("formed_overlap_counts") or []))}</code></div>'
        f'<div><b>Installed bond later broken</b><code>{html.escape(str(evidence.get("transient_bond_counts") or []))}</code></div>'
        f'<div><b>Replacement bonds</b><code>{html.escape(str(evidence.get("replacement_bond_counts") or []))}</code></div>'
        f'<div><b>Minimum graph distance</b><code>{html.escape(str(distances))}</code></div>'
        f'<div><b>Lineage invariant</b><code>{html.escape(str(bool(evidence.get("ambiguity_invariant"))).lower())}</code></div>'
        f'<div><b>Lineage status</b><code>{html.escape(str(item.get("lineage_status") or ""))}</code></div>'
        "</div>"
        f'<ul class="rationale">{rationale}</ul>'
        + (f'<ul class="warnings">{warnings}</ul>' if warnings else "")
    )


def _review_controls() -> str:
    return (
        '<div class="review">'
        '<label>Coupling classification<select class="review-coupling">'
        '<option value="unreviewed">Unreviewed</option>'
        '<option value="correct">Correct</option>'
        '<option value="too_strict">Too strict</option>'
        '<option value="too_permissive">Too permissive</option>'
        '<option value="uncertain">Uncertain</option></select></label>'
        '<label>Useful as one strategy?<select class="review-useful">'
        '<option value="unreviewed">Unreviewed</option>'
        '<option value="yes">Yes</option><option value="maybe">Maybe</option>'
        '<option value="no">No</option></select></label>'
        '<label>Suggested dependency<select class="review-dependency">'
        '<option value="unreviewed">Unreviewed</option>'
        '<option value="created_handle_consumed">Created handle consumed</option>'
        '<option value="activation_then_conversion">Activation then conversion</option>'
        '<option value="temporary_group_removed">Temporary group removed</option>'
        '<option value="continued_site_transformation">Continued site transformation</option>'
        '<option value="shared_local_environment">Shared local environment</option>'
        '<option value="independent_sites">Independent sites</option>'
        '<option value="lineage_ambiguous">Lineage ambiguous</option>'
        '<option value="unresolved">Unresolved</option></select></label>'
        '<label class="notes">Notes<textarea class="review-notes" rows="3" '
        'placeholder="Why these operations do or do not form one synthetic strategy…"></textarea></label>'
        "</div>"
    )


def _strategy_card(item: Mapping[str, Any], index: int) -> str:
    occurrence_id = str(item.get("occurrence_id") or "")
    relationship = str(item.get("relationship_class") or "unresolved")
    dependency = str(item.get("dependency_class") or relationship)
    admission = str(item.get("admission_class") or "review")
    route_id = str(item.get("source_route_id") or "")
    patent_id = str(item.get("patent_id") or "")
    first = str(item.get("first_reaction_smiles") or "")
    second = str(item.get("second_reaction_smiles") or "")
    intermediate = str(item.get("intermediate_smiles") or "")
    overall = str(item.get("overall_reaction_smiles") or "")
    search = " ".join(
        (
            occurrence_id,
            route_id,
            patent_id,
            relationship,
            dependency,
            admission,
            intermediate,
            overall,
            " ".join(str(value) for value in item.get("combined_edit_tokens") or ()),
        )
    ).lower()
    return (
        f'<article class="strategy-card admission-{html.escape(admission)}" '
        f'data-review-id="{html.escape(occurrence_id)}" '
        f'data-admission="{html.escape(admission)}" '
        f'data-relationship="{html.escape(relationship)}" '
        f'data-dependency="{html.escape(dependency)}" '
        f'data-route="{html.escape(route_id)}" data-search="{html.escape(search)}">'
        '<div class="card-head">'
        f'<span class="number">{index:02d}</span>'
        f'<span class="admission">{html.escape(admission)}</span>'
        f'<span class="relationship baseline">{html.escape(relationship.replace("_", " "))}</span>'
        '<span class="arrow">→</span>'
        f'<span class="dependency">{html.escape(dependency.replace("_", " "))}</span>'
        f'<b>score {float(item.get("coupling_score") or 0):.2f}</b>'
        f'<span class="source">{html.escape(patent_id)} · {html.escape(route_id)}</span>'
        "</div>"
        + (
            '<details class="overall" open><summary>Logical overall transformation</summary>'
            f'<div class="overall-reaction">{reaction_svg(overall)}</div>'
            f'<code>{html.escape(overall)}</code></details>'
            if overall
            else ""
        )
        + '<div class="sequence">'
        f'<div class="reaction"><h3>Physical step 1</h3>{reaction_svg(first)}'
        f'<ul>{_text_list(item.get("first_edit_tokens") or ())}</ul></div>'
        '<div class="intermediate"><h3>Carried intermediate</h3>'
        f'{molecule_svg(intermediate, width=360, height=230)}'
        f'<code>{html.escape(intermediate)}</code></div>'
        f'<div class="reaction"><h3>Physical step 2</h3>{reaction_svg(second)}'
        f'<ul>{_text_list(item.get("second_edit_tokens") or ())}</ul></div>'
        "</div>"
        '<details class="evidence" open><summary>Structural coupling evidence</summary>'
        + _evidence_panel(item)
        + '<div class="identity"><b>Typed strategy</b> '
        f'<code>{html.escape(str(item.get("typed_strategy_id") or ""))}</code>'
        f'<b>Occurrence</b> <code>{html.escape(occurrence_id)}</code></div></details>'
        + _review_controls()
        + "</article>"
    )


def _stylesheet() -> str:
    return """
:root{--ink:#17201d;--muted:#63716b;--line:#d4ddd8;--paper:#f2f5f3;--green:#176b4a;--amber:#9a6514;--red:#9e3e35;--blue:#346a8a}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}header{background:#193f34;color:white;padding:24px max(18px,calc((100% - 1560px)/2))}header h1{margin:0 0 5px;font-size:28px}header p{margin:0;max-width:1050px;color:#dceae4}.summary{display:flex;gap:8px;flex-wrap:wrap;margin-top:14px}.summary span{background:#2b5549;padding:7px 11px;border-radius:7px}.summary b{font-size:18px;margin-left:4px}main{max-width:1560px;margin:auto;padding:14px}.top-strategies,.comparison{background:white;border:1px solid var(--line);border-radius:8px;margin-bottom:10px;padding:9px}.top-strategies summary,.comparison summary{cursor:pointer;font-weight:650}.comparison p{color:var(--muted)}.top-strategies table,.comparison table{border-collapse:collapse;width:100%;margin-top:8px}.top-strategies th,.top-strategies td,.comparison th,.comparison td{text-align:left;vertical-align:top;border-top:1px solid var(--line);padding:6px}.comparison .arrow,.card-head .arrow{color:var(--blue);font-weight:700}.toolbar{position:sticky;top:0;z-index:10;display:flex;gap:9px;align-items:end;flex-wrap:wrap;background:white;border:1px solid var(--line);padding:10px;border-radius:9px;box-shadow:0 4px 15px #173c3020}.toolbar label,.review label{display:grid;gap:3px;color:var(--muted);font-size:12px}.toolbar input,.toolbar select,.review select,.review textarea{border:1px solid #b8c6bf;border-radius:6px;background:white;padding:7px}.toolbar input{width:280px}.toolbar button{border:0;border-radius:6px;background:#315e50;color:white;padding:8px 12px;cursor:pointer}.visible{margin-left:auto;color:var(--muted)}.strategy-card{background:white;border:1px solid var(--line);border-left:6px solid var(--amber);border-radius:10px;margin:15px 0;padding:12px;box-shadow:0 3px 12px #1f372e0c}.strategy-card.admission-strict{border-left-color:var(--green)}.strategy-card.admission-rejected{border-left-color:var(--red)}.strategy-card[data-useful="yes"]{box-shadow:inset 0 0 0 2px var(--green)}.strategy-card[data-useful="no"]{box-shadow:inset 0 0 0 2px var(--red)}.card-head{display:flex;gap:8px;align-items:center;flex-wrap:wrap}.number,.source{color:var(--muted)}.source{margin-left:auto}.admission,.relationship,.dependency{border-radius:999px;padding:4px 9px;background:#e6eeea}.relationship.baseline{background:#eef1ef;color:var(--muted)}.dependency{background:#dceee7;color:#124e37;font-weight:650}.sequence{display:grid;grid-template-columns:1fr 380px 1fr;gap:9px;margin:10px 0}.sequence>div,.overall{border:1px solid var(--line);border-radius:7px;padding:8px;overflow:auto}.sequence h3{margin:0 0 5px;font-size:14px}.reaction svg,.overall-reaction svg{width:100%;min-width:500px;height:auto;max-height:220px}.intermediate svg{width:100%;height:230px}.intermediate code,.overall code{display:block;overflow-wrap:anywhere}.reaction ul{margin:4px 0;padding-left:20px}.evidence{border:1px solid var(--line);border-radius:7px;padding:8px}.evidence summary,.overall summary{font-weight:650;cursor:pointer}.evidence-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:7px;margin-top:8px}.evidence-grid>div{background:#f6f8f7;border-radius:5px;padding:7px}.evidence-grid b,.evidence-grid code{display:block}.rationale{color:var(--green)}.warnings{color:var(--amber)}.identity{display:grid;grid-template-columns:110px 1fr;gap:4px 8px}.identity code{overflow-wrap:anywhere}.review{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;border-top:1px solid var(--line);padding-top:10px;margin-top:10px}.review .notes{grid-column:1/-1}.review textarea{width:100%;resize:vertical}@media(max-width:1100px){.sequence{grid-template-columns:1fr}.evidence-grid{grid-template-columns:repeat(2,1fr)}.toolbar{position:static}.visible{width:100%;margin:0}}@media print{.toolbar,.review{display:none}.strategy-card{break-inside:avoid}}
"""


def _script(storage_key: str, metadata: Mapping[str, Any]) -> str:
    return f"""
const key={json.dumps(storage_key)};let saved={{}};try{{saved=JSON.parse(localStorage.getItem(key)||'{{}}')}}catch(_){{saved={{}}}}const cards=[...document.querySelectorAll('.strategy-card')];for(const card of cards){{const id=card.dataset.reviewId,state=saved[id]||{{}};const controls={{coupling:'.review-coupling',useful:'.review-useful',dependency:'.review-dependency',notes:'.review-notes'}};for(const [name,selector] of Object.entries(controls)){{const el=card.querySelector(selector);el.value=state[name]||(name==='notes'?'':'unreviewed');el.addEventListener(name==='notes'?'input':'change',()=>{{const value={{}};for(const [field,sel] of Object.entries(controls))value[field]=card.querySelector(sel).value;saved[id]=value;card.dataset.useful=value.useful;try{{localStorage.setItem(key,JSON.stringify(saved))}}catch(_){{}}filter();}})}}card.dataset.useful=state.useful||'unreviewed';}}function filter(){{const query=document.querySelector('#search').value.toLowerCase(),admission=document.querySelector('#admission-filter').value,dependency=document.querySelector('#dependency-filter').value,useful=document.querySelector('#useful-filter').value;let visible=0;for(const card of cards){{const show=(!query||card.dataset.search.includes(query))&&(admission==='all'||card.dataset.admission===admission)&&(dependency==='all'||card.dataset.dependency===dependency)&&(useful==='all'||card.dataset.useful===useful);card.hidden=!show;if(show)visible++}}document.querySelector('#visible').textContent=`${{visible}} / ${{cards.length}} pairs`;}}for(const id of ['search','admission-filter','dependency-filter','useful-filter'])document.querySelector('#'+id).addEventListener('input',filter);document.querySelector('#export').addEventListener('click',()=>{{const reviews=cards.map(card=>({{occurrence_id:card.dataset.reviewId,admission_class:card.dataset.admission,relationship_class:card.dataset.relationship,dependency_class:card.dataset.dependency,...(saved[card.dataset.reviewId]||{{coupling:'unreviewed',useful:'unreviewed',dependency:'unreviewed',notes:''}})}}));const payload={{review_schema_version:'1.1',artifact:{json.dumps(dict(metadata), sort_keys=True)},reviews}};const blob=new Blob([JSON.stringify(payload,null,2)],{{type:'application/json'}}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='coupled-route-strategy-review.v2.json';a.click();URL.revokeObjectURL(a.href)}});filter();
"""


def render_coupled_route_strategy_html(
    report: Mapping[str, Any],
    *,
    title: str = "Coupled two-step strategy review",
) -> str:
    """Render the balanced mining sample as standalone review HTML."""

    sample = tuple(report.get("sample") or ())
    if not sample:
        raise ValueError("coupled route-strategy report has no review sample")
    counts = report.get("admission_counts") or {}
    cards = "".join(_strategy_card(item, index) for index, item in enumerate(sample, 1))
    dependencies = sorted(
        {str(item.get("dependency_class") or "") for item in sample}
    )
    dependency_options = "".join(
        f'<option value="{html.escape(value)}">{html.escape(value.replace("_", " ").title())}</option>'
        for value in dependencies
    )
    metadata = {
        "algorithm_version": report.get("algorithm_version"),
        "source_sha256": report.get("source_sha256"),
        "sample_seed": report.get("sample_seed"),
    }
    storage_key = "coupled-route-strategy-review:" + ":".join(
        (
            str(report.get("algorithm_version") or "unknown"),
            str(report.get("source_sha256") or "unknown"),
        )
    )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{html.escape(title)}</title><style>{_stylesheet()}</style></head><body>'
        f'<header><h1>{html.escape(title)}</h1>'
        '<p>Compare the v1 site-overlap relationship with the v2 causal dependency. Strict pairs now distinguish a created handle consumed later, activation followed by conversion, and continued transformation at one site. Physical steps and the intermediate remain explicit.</p>'
        '<div class="summary">'
        f'<span>Routes <b>{int(report.get("route_count") or 0)}</b></span>'
        f'<span>Pairs <b>{int(report.get("lineage_pair_count") or 0)}</b></span>'
        f'<span>Strict <b>{int(counts.get("strict") or 0)}</b></span>'
        f'<span>Review <b>{int(counts.get("review") or 0)}</b></span>'
        f'<span>Rejected <b>{int(counts.get("rejected") or 0)}</b></span>'
        f'<span>Review sample <b>{len(sample)}</b></span></div></header>'
        '<main>'
        + _comparison_table(report)
        + _strategy_summary_table(report)
        + '<section class="toolbar"><label>Search<input id="search" type="search" placeholder="Patent, route, intermediate, edit…"></label>'
        '<label>Admission<select id="admission-filter"><option value="all">All</option><option value="strict">Strict</option><option value="review">Review</option><option value="rejected">Rejected control</option></select></label>'
        '<label>Dependency<select id="dependency-filter"><option value="all">All</option>'
        + dependency_options
        + '</select></label><label>Useful review<select id="useful-filter"><option value="all">All</option><option value="unreviewed">Unreviewed</option><option value="yes">Yes</option><option value="maybe">Maybe</option><option value="no">No</option></select></label>'
        '<button id="export" type="button">Export review JSON</button><span class="visible" id="visible"></span></section>'
        + cards
        + f'</main><script>{_script(storage_key, metadata)}</script></body></html>'
    )


def write_coupled_route_strategy_html(
    source_report: Mapping[str, Any] | str | Path,
    output_html: str | Path,
    *,
    title: str = "Coupled two-step strategy review",
) -> dict[str, Any]:
    """Write the self-contained HTML review and return its summary."""

    report = (
        load_coupled_route_strategy_report(source_report)
        if isinstance(source_report, (str, Path))
        else dict(source_report)
    )
    document = render_coupled_route_strategy_html(report, title=title)
    destination = Path(output_html)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(document, encoding="utf-8", newline="\n")
    return {
        "review_version": COUPLED_ROUTE_STRATEGY_REVIEW_VERSION,
        "sample_count": len(report.get("sample") or ()),
        "output_html": str(destination.resolve()),
        "html_bytes": destination.stat().st_size,
    }


__all__ = [
    "load_coupled_route_strategy_report",
    "render_coupled_route_strategy_html",
    "write_coupled_route_strategy_html",
]
