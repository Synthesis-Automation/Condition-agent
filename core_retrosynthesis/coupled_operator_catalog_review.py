"""Self-contained graphical catalog for executable v1 operator pairs."""

from __future__ import annotations

import html
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .coupled_strategy_evaluation import (
    CoupledStrategyEvaluationCase,
    PromotedV1OperatorPair,
)
from .generic_models import GenericGraphOperator, GenericTemplateLibrary
from .html_report import reaction_svg


COUPLED_OPERATOR_CATALOG_REVIEW_VERSION = "v1_coupled_operator_catalog_review.v2"


def _display_name(value: str) -> str:
    return value.replace("_", " ").strip().title()


def _compact_case(
    cases: Iterable[CoupledStrategyEvaluationCase],
) -> CoupledStrategyEvaluationCase | None:
    values = tuple(cases)
    if not values:
        return None
    return min(
        values,
        key=lambda item: (
            len(item.observed_first_reaction_smiles)
            + len(item.observed_second_reaction_smiles),
            item.patent_id,
            item.case_id,
        ),
    )


def _operator_summary(
    operator: GenericGraphOperator | None,
    *,
    label: str,
) -> str:
    if operator is None:
        return (
            f'<div class="operator-summary missing"><strong>{label}</strong>'
            "<span>Operator unavailable</span></div>"
        )
    levels = " / ".join(sorted(set(operator.abstraction_levels))) or "—"
    edits = "".join(
        f'<span class="edit-token">{html.escape(token)}</span>'
        for token in operator.edit_tokens[:6]
    )
    if len(operator.edit_tokens) > 6:
        edits += (
            f'<span class="edit-token more">+{len(operator.edit_tokens) - 6}'
            " more</span>"
        )
    return (
        '<div class="operator-summary">'
        f"<strong>{html.escape(label)}</strong>"
        f"<code>{html.escape(operator.operator_id)}</code>"
        '<div class="operator-stats">'
        f"<span>{levels}</span>"
        f"<span>{operator.observation_support} observations</span>"
        f"<span>{operator.independent_reference_support} references</span>"
        "</div>"
        f'<div class="edit-list">{edits}</div>'
        "</div>"
    )


def _step_graphic(
    reaction_smiles: str,
    *,
    number: int,
    caption: str,
) -> str:
    return (
        '<section class="physical-step">'
        '<div class="step-heading">'
        f"<span>Physical step {number}</span><strong>{html.escape(caption)}</strong>"
        "</div>"
        f'<div class="reaction-graphic">{reaction_svg(reaction_smiles, sub_image_size=(500, 300))}</div>'
        "</section>"
    )


def _pair_card(
    strategy: PromotedV1OperatorPair,
    cases: Sequence[CoupledStrategyEvaluationCase],
    operators: Mapping[str, GenericGraphOperator],
    index: int,
) -> str:
    example = _compact_case(cases)
    first_operator = operators.get(strategy.first_operator_id)
    second_operator = operators.get(strategy.second_operator_id)
    dependency_counts = dict(strategy.v2_dependency_counts)
    dependency_text = " · ".join(
        f"{_display_name(name)} ({count})"
        for name, count in sorted(dependency_counts.items())
    )
    search_value = " ".join(
        (
            strategy.strategy_id,
            strategy.relationship_class,
            strategy.first_operator_id,
            strategy.second_operator_id,
            *dependency_counts,
            *(first_operator.edit_tokens if first_operator else ()),
            *(second_operator.edit_tokens if second_operator else ()),
        )
    ).casefold()
    if example is None:
        route_graphic = (
            '<div class="render-error">No held-out observed example available</div>'
        )
        intermediate = "—"
        terminal = "—"
        target = "—"
        patent_id = "—"
        case_id = "—"
        first_reaction = ""
        second_reaction = ""
    else:
        intermediate = example.expected_intermediate_smiles
        terminal = example.expected_terminal_precursor_smiles
        target = example.target_smiles
        patent_id = example.patent_id
        case_id = example.case_id
        first_reaction = example.observed_first_reaction_smiles
        second_reaction = example.observed_second_reaction_smiles
        route_graphic = (
            '<div class="two-step-route">'
            + _step_graphic(
                first_reaction,
                number=1,
                caption="Terminal precursors → intermediate",
            )
            + '<div class="route-bridge"><span>then</span>'
            f"<code>{html.escape(intermediate)}</code></div>"
            + _step_graphic(
                second_reaction,
                number=2,
                caption="Intermediate → final product",
            )
            + "</div>"
        )
    return (
        '<article class="pair-card" '
        f'data-search="{html.escape(search_value, quote=True)}" '
        f'data-relationship="{html.escape(strategy.relationship_class, quote=True)}" '
        f'data-dependencies="{html.escape(" ".join(dependency_counts), quote=True)}">'
        '<header class="pair-heading">'
        f'<span class="pair-number">{index:03d}</span>'
        '<div><span class="eyebrow">EXECUTABLE TWO-STEP OPERATOR PAIR</span>'
        f"<h2>{html.escape(strategy.strategy_id)}</h2>"
        f"<p>{html.escape(_display_name(strategy.relationship_class))}</p></div>"
        '<div class="support-badge">'
        f"<strong>{len(strategy.training_patent_ids)}</strong><span>training patents</span>"
        f"<strong>{strategy.training_occurrence_count}</strong><span>occurrences</span>"
        "</div>"
        "</header>"
        f'<div class="dependency-line"><strong>V2 review labels</strong><span>{html.escape(dependency_text or "None")}</span></div>'
        f"{route_graphic}"
        '<div class="operator-pair">'
        f"{_operator_summary(first_operator, label='Step 1 operator')}"
        '<div class="pair-arrow">→</div>'
        f"{_operator_summary(second_operator, label='Step 2 operator')}"
        "</div>"
        "<details><summary>Observed example and full identities</summary><dl>"
        f"<div><dt>Terminal precursors</dt><dd><code>{html.escape(terminal)}</code></dd></div>"
        f"<div><dt>Intermediate</dt><dd><code>{html.escape(intermediate)}</code></dd></div>"
        f"<div><dt>Final product</dt><dd><code>{html.escape(target)}</code></dd></div>"
        f"<div><dt>First reaction</dt><dd><code>{html.escape(first_reaction)}</code></dd></div>"
        f"<div><dt>Second reaction</dt><dd><code>{html.escape(second_reaction)}</code></dd></div>"
        f"<div><dt>Example patent</dt><dd><code>{html.escape(patent_id)}</code></dd></div>"
        f"<div><dt>Case ID</dt><dd><code>{html.escape(case_id)}</code></dd></div>"
        f"<div><dt>Training patents</dt><dd><code>{html.escape(' · '.join(strategy.training_patent_ids))}</code></dd></div>"
        "</dl></details>"
        "</article>"
    )


def render_v1_operator_pair_catalog_html(
    strategies: Sequence[PromotedV1OperatorPair],
    cases: Sequence[CoupledStrategyEvaluationCase],
    library: GenericTemplateLibrary,
    *,
    title: str = "Executable two-step operator-pair catalog",
    route_core_source: str | Path | None = None,
    library_source: str | Path | None = None,
    capability_gap_count: int = 0,
) -> str:
    """Render every supplied v1 pair with two connected reaction graphics."""

    ordered = tuple(sorted(strategies, key=lambda item: item.strategy_id))
    if len({item.strategy_id for item in ordered}) != len(ordered):
        raise ValueError("two-step operator catalog requires unique strategy IDs")
    cases_by_strategy: dict[str, list[CoupledStrategyEvaluationCase]] = defaultdict(
        list
    )
    for case in cases:
        cases_by_strategy[case.strategy_id].append(case)
    operators = {item.operator_id: item for item in library.operators}
    cards = "".join(
        _pair_card(
            strategy,
            tuple(cases_by_strategy.get(strategy.strategy_id, ())),
            operators,
            index,
        )
        for index, strategy in enumerate(ordered, 1)
    )
    relationship_counts = Counter(item.relationship_class for item in ordered)
    dependency_values = sorted(
        {name for strategy in ordered for name, _count in strategy.v2_dependency_counts}
    )
    dependency_options = "".join(
        f'<option value="{html.escape(value, quote=True)}">{html.escape(_display_name(value))}</option>'
        for value in dependency_values
    )
    metadata = {
        "review_version": COUPLED_OPERATOR_CATALOG_REVIEW_VERSION,
        "rendered_pair_count": len(ordered),
        "relationship_counts": dict(sorted(relationship_counts.items())),
        "capability_gap_count": capability_gap_count,
        "ordering": "strategy_id_ascending",
        "graphic_source": "shortest_patent_disjoint_heldout_example_per_pair",
        "graphic_layout": "one_physical_reaction_per_row",
        "minimum_training_patents": 2,
    }
    route_name = (
        Path(route_core_source).name if route_core_source else "in-memory routes"
    )
    library_name = Path(library_source).name if library_source else "in-memory library"
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{ font-family: Inter, ui-sans-serif, system-ui, sans-serif; color: #193126; background: #eef3ef; }}
* {{ box-sizing: border-box; }} body {{ margin: 0; }}
.shell {{ width: min(1580px,100%); margin: 0 auto; padding: 24px; }}
.page-heading {{ display: flex; justify-content: space-between; gap: 24px; align-items: flex-start; margin-bottom: 18px; }}
h1 {{ margin: 0; color: #153829; font: 500 2rem Georgia,serif; }}
.page-heading p {{ max-width: 930px; margin: 8px 0 0; color: #5d7166; line-height: 1.5; }}
.summary {{ display: grid; grid-template-columns: repeat(2,minmax(105px,1fr)); gap: 1px; overflow: hidden; border: 1px solid #cedbd3; border-radius: 10px; background: #cedbd3; }}
.summary span {{ padding: 10px 14px; background: #fff; text-align: center; }}
.summary strong {{ display: block; color: #28533e; font-size: 1.05rem; }}
.summary small {{ color: #718078; font-size: .62rem; text-transform: uppercase; }}
.notice {{ margin-bottom: 15px; padding: 11px 13px; border: 1px solid #e2d2aa; border-radius: 9px; color: #6f541e; background: #fff8e8; font-size: .75rem; line-height: 1.45; }}
.toolbar {{ position: sticky; z-index: 5; top: 0; display: grid; grid-template-columns: minmax(260px,1fr) 220px 250px auto; gap: 10px; margin-bottom: 18px; padding: 12px; border: 1px solid #cedbd3; border-radius: 11px; background: rgba(255,255,255,.97); box-shadow: 0 8px 24px rgba(32,62,47,.08); }}
input,select {{ min-height: 40px; border: 1px solid #c6d5cc; border-radius: 8px; padding: 0 11px; color: #244536; background: #fff; font: inherit; }}
.visible-count {{ align-self: center; color: #61756a; font-size: .76rem; text-align: right; }}
.catalog {{ display: grid; gap: 17px; }}
.pair-card {{ overflow: hidden; border: 1px solid #cfdbd3; border-radius: 12px; background: #fff; box-shadow: 0 8px 25px rgba(37,64,51,.06); }}
.pair-card[hidden] {{ display: none; }}
.pair-heading {{ display: grid; grid-template-columns: 50px minmax(0,1fr) minmax(220px,auto); gap: 12px; align-items: start; padding: 15px 17px 10px; }}
.pair-number {{ display: grid; width: 48px; height: 40px; place-items: center; border-radius: 8px; color: #fff; background: #2c6d50; font: 800 .78rem Consolas,monospace; }}
.eyebrow {{ color: #54806b; font-size: .6rem; font-weight: 850; letter-spacing: .13em; }}
h2 {{ overflow-wrap: anywhere; margin: 3px 0 4px; color: #264d3a; font: 700 .76rem Consolas,monospace; }}
.pair-heading p {{ margin: 0; color: #687b70; font-size: .7rem; }}
.support-badge {{ display: grid; grid-template-columns: auto 1fr; gap: 2px 7px; padding: 7px 10px; border: 1px solid #d5e1d9; border-radius: 8px; background: #f7faf8; }}
.support-badge strong {{ color: #28533e; font-size: .78rem; text-align: right; }} .support-badge span {{ color: #728179; font-size: .64rem; }}
.dependency-line {{ display: flex; gap: 10px; padding: 8px 17px; border-block: 1px solid #e1e8e3; color: #687a70; background: #f8faf8; font-size: .67rem; }}
.dependency-line strong {{ color: #456653; }}
.two-step-route {{ display: grid; grid-template-columns: minmax(0,1fr); gap: 10px; align-items: stretch; padding: 15px 17px 12px; }}
.physical-step {{ min-width: 0; overflow: hidden; border: 1px solid #d6e0d9; border-radius: 9px; background: #fff; }}
.step-heading {{ display: flex; justify-content: space-between; gap: 12px; padding: 9px 12px; border-bottom: 1px solid #e0e7e2; color: #687a70; background: #f7faf8; font-size: .72rem; }}
.step-heading strong {{ color: #3f614f; }} .reaction-graphic svg {{ display: block; width: 100%; height: clamp(240px,22vw,340px); }}
.route-bridge {{ display: flex; min-width: 0; min-height: 58px; flex-direction: column; justify-content: center; gap: 6px; padding: 4px 18px; color: #687a70; text-align: center; }}
.route-bridge span {{ align-self: center; padding: 4px 10px; border-radius: 99px; color: #fff; background: #5d826e; font-size: .63rem; font-weight: 800; text-transform: uppercase; }}
.route-bridge code {{ overflow-wrap: anywhere; color: #3d5d4d; font-size: .61rem; }}
.operator-pair {{ display: grid; grid-template-columns: minmax(0,1fr) 34px minmax(0,1fr); gap: 8px; align-items: center; padding: 0 17px 14px; }}
.operator-summary {{ min-width: 0; padding: 10px; border: 1px solid #d5dfd9; border-radius: 8px; background: #f8faf8; }}
.operator-summary > strong {{ display: block; margin-bottom: 4px; color: #3e604e; font-size: .67rem; }}
.operator-summary > code {{ display: block; overflow-wrap: anywhere; color: #355743; font-size: .61rem; }}
.operator-stats {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 6px; color: #718078; font-size: .61rem; }}
.pair-arrow {{ color: #4f7863; font-size: 1.3rem; text-align: center; }}
.edit-list {{ display: flex; flex-wrap: wrap; gap: 4px; margin-top: 7px; }}
.edit-token {{ padding: 3px 6px; border: 1px solid #e3d6b9; border-radius: 5px; color: #695121; background: #fff8e9; font: .58rem Consolas,monospace; }}
.edit-token.more {{ color: #68786f; border-color: #d5dfd9; background: #f6f8f6; }}
.render-error {{ display: grid; min-height: 180px; place-items: center; color: #98533f; }}
details {{ border-top: 1px solid #e0e7e2; }} summary {{ padding: 10px 17px; color: #456653; cursor: pointer; font-size: .7rem; font-weight: 750; }}
dl {{ display: grid; gap: 7px; margin: 0; padding: 0 17px 15px; }} dl div {{ display: grid; grid-template-columns: 145px minmax(0,1fr); gap: 9px; }}
dt {{ color: #748279; font-size: .64rem; }} dd {{ min-width: 0; margin: 0; color: #425f50; font-size: .65rem; }} code {{ overflow-wrap: anywhere; font-family: Consolas,monospace; }}
.metadata {{ margin-top: 18px; padding: 12px; border: 1px solid #d3ded7; border-radius: 9px; color: #66786e; background: #f8faf8; font-size: .67rem; white-space: pre-wrap; }}
@media (max-width:1050px) {{ .operator-pair {{ grid-template-columns: 1fr; }} .pair-arrow {{ transform: rotate(90deg); }} .toolbar {{ grid-template-columns: 1fr 1fr; }} }}
@media (max-width:650px) {{ .shell {{ padding: 12px; }} .page-heading {{ flex-direction: column; }} .toolbar,.pair-heading {{ grid-template-columns: 1fr; }} .visible-count {{ text-align: left; }} }}
@media print {{ .toolbar {{ display:none; }} .pair-card {{ break-inside:avoid; }} }}
</style>
</head>
<body><main class="shell">
<header class="page-heading"><div><span class="eyebrow">ROUTE-DERIVED V1 STRATEGIES</span><h1>{html.escape(title)}</h1><p>All {len(ordered)} recurrent, patent-disjoint v1 operator pairs for which both physical one-step operators are executable in the validated-departures route-only library. Each card uses the shortest held-out observed route example for legibility.</p></div><div class="summary"><span><strong>{len(ordered)}</strong><small>executable pairs</small></span><span><strong>{capability_gap_count}</strong><small>coverage gaps</small></span><span><strong>{relationship_counts.get("same_site_coupled", 0)}</strong><small>same-site</small></span><span><strong>{relationship_counts.get("handle_progression", 0)}</strong><small>handle progression</small></span></div></header>
<div class="notice"><strong>Experimental review catalog.</strong> Both physical operators are executable and structurally supported, but the pair is not automatically a production-approved synthetic strategy. Review selectivity, intermediate stability, and condition compatibility.</div>
<section class="toolbar"><input id="search" type="search" placeholder="Search pair ID, operator ID, edit, or dependency"><select id="relationship"><option value="">All relationships</option><option value="handle_progression">Handle progression</option><option value="same_site_coupled">Same-site coupled</option></select><select id="dependency"><option value="">All V2 review labels</option>{dependency_options}</select><span id="visible-count" class="visible-count">{len(ordered)} pairs visible</span></section>
<section class="catalog">{cards}</section>
<pre class="metadata">Route source: {html.escape(route_name)}
Operator library: {html.escape(library_name)}
{html.escape(json.dumps(metadata, indent=2, sort_keys=True))}</pre>
</main><script>
const search=document.getElementById('search'); const relationship=document.getElementById('relationship'); const dependency=document.getElementById('dependency'); const count=document.getElementById('visible-count'); const cards=[...document.querySelectorAll('.pair-card')];
function filterCards() {{ const query=search.value.trim().toLowerCase(); const selectedRelationship=relationship.value; const selectedDependency=dependency.value; let visible=0; for (const card of cards) {{ const matchesText=!query||card.dataset.search.includes(query); const matchesRelationship=!selectedRelationship||card.dataset.relationship===selectedRelationship; const matchesDependency=!selectedDependency||card.dataset.dependencies.split(' ').includes(selectedDependency); card.hidden=!(matchesText&&matchesRelationship&&matchesDependency); if(!card.hidden) visible+=1; }} count.textContent=`${{visible}} pair${{visible===1?'':'s'}} visible`; }}
search.addEventListener('input',filterCards); relationship.addEventListener('change',filterCards); dependency.addEventListener('change',filterCards);
</script></body></html>"""


def write_v1_operator_pair_catalog_html(
    strategies: Sequence[PromotedV1OperatorPair],
    cases: Sequence[CoupledStrategyEvaluationCase],
    library: GenericTemplateLibrary,
    output_path: str | Path,
    *,
    title: str = "Executable two-step operator-pair catalog",
    route_core_source: str | Path | None = None,
    library_source: str | Path | None = None,
    capability_gap_count: int = 0,
) -> dict[str, object]:
    """Write the complete self-contained graphical pair catalog."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    page = render_v1_operator_pair_catalog_html(
        strategies,
        cases,
        library,
        title=title,
        route_core_source=route_core_source,
        library_source=library_source,
        capability_gap_count=capability_gap_count,
    )
    destination.write_text(page, encoding="utf-8", newline="\n")
    return {
        "output_html": str(destination.resolve()),
        "html_bytes": destination.stat().st_size,
        "rendered_pair_count": len(strategies),
        "rendered_reaction_graphic_count": 2 * len(strategies),
        "capability_gap_count": capability_gap_count,
        "review_version": COUPLED_OPERATOR_CATALOG_REVIEW_VERSION,
    }


__all__ = [
    "COUPLED_OPERATOR_CATALOG_REVIEW_VERSION",
    "render_v1_operator_pair_catalog_html",
    "write_v1_operator_pair_catalog_html",
]
