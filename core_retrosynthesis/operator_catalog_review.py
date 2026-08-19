"""Self-contained graphical review catalog for generic graph operators."""

from __future__ import annotations

import html
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

from .generic_models import (
    GenericCoreTemplate,
    GenericGraphOperator,
    GenericTemplateLibrary,
    GenericTemplatePrecedent,
)
from .html_report import reaction_svg


OPERATOR_CATALOG_REVIEW_VERSION = "generic_operator_catalog_review.v1"


def _display_name(value: str) -> str:
    return value.replace("_", " ").strip().title()


def _compact_precedent(
    templates: Iterable[GenericCoreTemplate],
) -> tuple[GenericCoreTemplate, GenericTemplatePrecedent] | None:
    """Return a deterministic compact observed example for graphic review."""

    candidates = [
        (template, precedent)
        for template in templates
        for precedent in template.precedents
    ]
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda item: (
            len(item[1].precursor_smiles) + len(item[1].product_smiles),
            -item[0].independent_reference_support,
            -item[0].observation_support,
            item[0].template_id,
            item[1].reaction_id,
        ),
    )


def _edit_badges(edit_tokens: Sequence[str]) -> str:
    values = "".join(
        f'<span class="edit-token">{html.escape(token)}</span>'
        for token in edit_tokens[:8]
    )
    if len(edit_tokens) > 8:
        values += f'<span class="edit-token more">+{len(edit_tokens) - 8} more</span>'
    return values or '<span class="edit-token none">No edit token</span>'


def _operator_card(
    operator: GenericGraphOperator,
    templates: Sequence[GenericCoreTemplate],
    completion_group_count: int,
    index: int,
) -> str:
    example = _compact_precedent(templates)
    transformation_kinds = tuple(
        sorted(
            {
                template.transformation_kind
                for template in templates
                if template.transformation_kind
            }
        )
    )
    levels = tuple(sorted(set(operator.abstraction_levels)))
    search_value = " ".join(
        (
            operator.operator_id,
            *operator.edit_tokens,
            *operator.named_annotations,
            *transformation_kinds,
            *levels,
        )
    ).casefold()
    if example is None:
        graphic = '<div class="render-error">No observed precedent available</div>'
        reaction_smiles = ""
        reaction_id = "—"
        reference_id = "—"
        template_id = "—"
    else:
        template, precedent = example
        reaction_smiles = f"{precedent.precursor_smiles}>>{precedent.product_smiles}"
        graphic = reaction_svg(
            reaction_smiles,
            sub_image_size=(220, 145),
        )
        reaction_id = precedent.reaction_id
        reference_id = precedent.reference_id
        template_id = template.template_id
    kinds = (
        " · ".join(_display_name(value) for value in transformation_kinds)
        if transformation_kinds
        else "Unnamed structural transformation"
    )
    annotation_badges = "".join(
        f'<span class="annotation">{html.escape(value)}</span>'
        for value in operator.named_annotations
    )
    return (
        '<article class="operator-card" '
        f'data-search="{html.escape(search_value, quote=True)}" '
        f'data-levels="{html.escape(" ".join(levels), quote=True)}">'
        '<header class="operator-heading">'
        f'<span class="operator-number">{index:03d}</span>'
        '<div><span class="eyebrow">GENERIC GRAPH OPERATOR</span>'
        f"<h2>{html.escape(operator.operator_id)}</h2>"
        f"<p>{html.escape(kinds)}</p></div>"
        "</header>"
        '<div class="operator-metrics">'
        f"<span><strong>{operator.observation_support}</strong> observations</span>"
        f"<span><strong>{operator.independent_reference_support}</strong> references</span>"
        f"<span><strong>{len(templates)}</strong> templates</span>"
        f"<span><strong>{len(operator.realization_ids)}</strong> realizations</span>"
        f"<span><strong>{completion_group_count}</strong> completion groups</span>"
        "</div>"
        f'<div class="level-row">{"".join(f"<span>{html.escape(level)}</span>" for level in levels)}</div>'
        '<section class="graphic">'
        '<div class="graphic-label"><strong>Compact observed example</strong>'
        "<span>forward direction</span></div>"
        f"{graphic}</section>"
        f'<div class="edit-list">{_edit_badges(operator.edit_tokens)}</div>'
        f'<div class="annotation-list">{annotation_badges}</div>'
        "<details><summary>Operator and precedent details</summary>"
        "<dl>"
        f"<div><dt>Example reaction</dt><dd><code>{html.escape(reaction_smiles)}</code></dd></div>"
        f"<div><dt>Reaction ID</dt><dd><code>{html.escape(reaction_id)}</code></dd></div>"
        f"<div><dt>Reference ID</dt><dd><code>{html.escape(reference_id)}</code></dd></div>"
        f"<div><dt>Template ID</dt><dd><code>{html.escape(template_id)}</code></dd></div>"
        f"<div><dt>Operator signature</dt><dd><code>{html.escape(operator.operator_signature)}</code></dd></div>"
        "</dl></details>"
        "</article>"
    )


def render_generic_operator_catalog_html(
    library: GenericTemplateLibrary,
    *,
    limit: int = 100,
    title: str = "Generic operator graphical catalog",
    library_source: str | Path | None = None,
) -> str:
    """Render a searchable catalog for the first deterministic operators."""

    if limit < 1:
        raise ValueError("operator catalog limit must be positive")
    templates_by_operator: dict[str, list[GenericCoreTemplate]] = defaultdict(list)
    for template in library.templates:
        templates_by_operator[template.operator_id].append(template)
    completion_counts: dict[str, int] = defaultdict(int)
    for group in library.completion_groups:
        completion_counts[group.operator_id] += 1
    selected = tuple(
        sorted(library.operators, key=lambda item: item.operator_id)[:limit]
    )
    cards = "".join(
        _operator_card(
            operator,
            tuple(templates_by_operator.get(operator.operator_id, ())),
            completion_counts.get(operator.operator_id, 0),
            index,
        )
        for index, operator in enumerate(selected, 1)
    )
    source_name = Path(library_source).name if library_source else "in-memory library"
    metadata = {
        "review_version": OPERATOR_CATALOG_REVIEW_VERSION,
        "library_schema_version": library.schema_version,
        "library_operator_count": len(library.operators),
        "rendered_operator_count": len(selected),
        "ordering": "operator_id_ascending",
        "graphic_source": "shortest_observed_precedent_per_operator",
    }
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{ font-family: Inter, ui-sans-serif, system-ui, sans-serif; color: #193126; background: #eef3ef; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; }}
.shell {{ width: min(1500px, 100%); margin: 0 auto; padding: 24px; }}
.page-heading {{ display: flex; justify-content: space-between; gap: 24px; align-items: flex-start; margin-bottom: 18px; }}
h1 {{ margin: 0; color: #153829; font-family: Georgia, serif; font-size: 2rem; font-weight: 500; }}
.page-heading p {{ max-width: 850px; margin: 8px 0 0; color: #5d7166; line-height: 1.5; }}
.summary {{ display: flex; gap: 1px; overflow: hidden; border: 1px solid #cedbd3; border-radius: 10px; background: #cedbd3; }}
.summary span {{ min-width: 110px; padding: 10px 14px; background: #fff; text-align: center; }}
.summary strong {{ display: block; color: #28533e; font-size: 1.05rem; }}
.summary small {{ color: #718078; font-size: .66rem; text-transform: uppercase; }}
.toolbar {{ position: sticky; z-index: 5; top: 0; display: grid; grid-template-columns: minmax(260px, 1fr) 180px auto; gap: 10px; margin-bottom: 18px; padding: 12px; border: 1px solid #cedbd3; border-radius: 11px; background: rgba(255,255,255,.96); box-shadow: 0 8px 24px rgba(32,62,47,.08); }}
input, select, button {{ min-height: 40px; border: 1px solid #c6d5cc; border-radius: 8px; padding: 0 11px; color: #244536; background: #fff; font: inherit; }}
button {{ cursor: pointer; font-weight: 700; }}
.visible-count {{ align-self: center; color: #61756a; font-size: .78rem; text-align: right; }}
.catalog {{ display: grid; grid-template-columns: repeat(2, minmax(0,1fr)); gap: 15px; }}
.operator-card {{ min-width: 0; overflow: hidden; border: 1px solid #cfdbd3; border-radius: 12px; background: #fff; box-shadow: 0 8px 25px rgba(37,64,51,.06); }}
.operator-card[hidden] {{ display: none; }}
.operator-heading {{ display: grid; grid-template-columns: 48px minmax(0,1fr); gap: 12px; align-items: start; padding: 15px 16px 10px; }}
.operator-number {{ display: grid; width: 46px; height: 38px; place-items: center; border-radius: 8px; color: #fff; background: #2c6d50; font-family: Consolas, monospace; font-weight: 800; }}
.eyebrow {{ color: #54806b; font-size: .6rem; font-weight: 850; letter-spacing: .13em; }}
h2 {{ overflow-wrap: anywhere; margin: 3px 0 4px; color: #264d3a; font: 700 .75rem Consolas, monospace; }}
.operator-heading p {{ margin: 0; color: #687b70; font-size: .7rem; }}
.operator-metrics {{ display: flex; flex-wrap: wrap; gap: 1px; margin: 0 16px 9px; overflow: hidden; border: 1px solid #dbe5de; border-radius: 8px; background: #dbe5de; }}
.operator-metrics span {{ flex: 1 1 105px; padding: 7px 8px; color: #708078; background: #f8faf8; font-size: .62rem; text-align: center; }}
.operator-metrics strong {{ color: #315a46; }}
.level-row {{ display: flex; gap: 5px; margin: 0 16px 10px; }}
.level-row span, .annotation {{ padding: 3px 7px; border-radius: 99px; color: #315a46; background: #e8f2ec; font-size: .61rem; font-weight: 800; }}
.graphic {{ margin: 0 16px 11px; overflow: hidden; border: 1px solid #d6e0d9; border-radius: 9px; background: #fff; }}
.graphic-label {{ display: flex; justify-content: space-between; padding: 7px 9px; color: #687a70; border-bottom: 1px solid #e0e7e2; background: #f7faf8; font-size: .64rem; }}
.graphic svg {{ display: block; width: 100%; height: auto; min-height: 145px; }}
.render-error {{ display: grid; min-height: 145px; place-items: center; color: #98533f; }}
.edit-list, .annotation-list {{ display: flex; flex-wrap: wrap; gap: 5px; margin: 0 16px 10px; }}
.edit-token {{ padding: 4px 7px; border: 1px solid #e3d6b9; border-radius: 6px; color: #695121; background: #fff8e9; font: .61rem Consolas, monospace; }}
.edit-token.more, .edit-token.none {{ color: #68786f; border-color: #d5dfd9; background: #f6f8f6; }}
details {{ border-top: 1px solid #e0e7e2; }}
summary {{ padding: 10px 16px; color: #456653; cursor: pointer; font-size: .7rem; font-weight: 750; }}
dl {{ display: grid; gap: 7px; margin: 0; padding: 0 16px 14px; }}
dl div {{ display: grid; grid-template-columns: 120px minmax(0,1fr); gap: 9px; }}
dt {{ color: #748279; font-size: .64rem; }} dd {{ min-width: 0; margin: 0; color: #425f50; font-size: .65rem; }}
code {{ overflow-wrap: anywhere; font-family: Consolas, monospace; }}
.metadata {{ margin-top: 18px; padding: 12px; border: 1px solid #d3ded7; border-radius: 9px; color: #66786e; background: #f8faf8; font-size: .67rem; white-space: pre-wrap; }}
@media (max-width: 950px) {{ .catalog {{ grid-template-columns: 1fr; }} .page-heading {{ flex-direction: column; }} }}
@media (max-width: 620px) {{ .shell {{ padding: 12px; }} .toolbar {{ grid-template-columns: 1fr; }} .visible-count {{ text-align: left; }} }}
@media print {{ .toolbar {{ display: none; }} .catalog {{ grid-template-columns: 1fr; }} .operator-card {{ break-inside: avoid; }} }}
</style>
</head>
<body>
<main class="shell">
<header class="page-heading"><div><span class="eyebrow">ROUTE-STEP OPERATOR LIBRARY</span><h1>{html.escape(title)}</h1><p>First {len(selected)} operators ordered by stable operator ID. Each graphic is the shortest real observed precedent attached to that operator; it illustrates one realization and is not the complete generic scope.</p></div><div class="summary"><span><strong>{len(selected)}</strong><small>shown</small></span><span><strong>{len(library.operators)}</strong><small>library total</small></span></div></header>
<section class="toolbar"><input id="search" type="search" placeholder="Search operator ID, edit, transformation, or annotation"><select id="level"><option value="">All abstraction levels</option><option value="L0">L0</option><option value="L1">L1</option><option value="L2">L2</option></select><span id="visible-count" class="visible-count">{len(selected)} operators visible</span></section>
<section id="catalog" class="catalog">{cards}</section>
<pre class="metadata">Source: {html.escape(source_name)}
{html.escape(json.dumps(metadata, indent=2, sort_keys=True))}</pre>
</main>
<script>
const search = document.getElementById('search');
const level = document.getElementById('level');
const count = document.getElementById('visible-count');
const cards = [...document.querySelectorAll('.operator-card')];
function filterCards() {{
  const query = search.value.trim().toLowerCase();
  const selectedLevel = level.value;
  let visible = 0;
  for (const card of cards) {{
    const matchesText = !query || card.dataset.search.includes(query);
    const matchesLevel = !selectedLevel || card.dataset.levels.split(' ').includes(selectedLevel);
    card.hidden = !(matchesText && matchesLevel);
    if (!card.hidden) visible += 1;
  }}
  count.textContent = `${{visible}} operator${{visible === 1 ? '' : 's'}} visible`;
}}
search.addEventListener('input', filterCards);
level.addEventListener('change', filterCards);
</script>
</body>
</html>
"""


def write_generic_operator_catalog_html(
    library: GenericTemplateLibrary,
    output_path: str | Path,
    *,
    limit: int = 100,
    title: str = "Generic operator graphical catalog",
    library_source: str | Path | None = None,
) -> dict[str, object]:
    """Write the self-contained graphical operator catalog."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    page = render_generic_operator_catalog_html(
        library,
        limit=limit,
        title=title,
        library_source=library_source,
    )
    destination.write_text(page, encoding="utf-8", newline="\n")
    rendered_count = min(limit, len(library.operators))
    return {
        "output_html": str(destination.resolve()),
        "html_bytes": destination.stat().st_size,
        "rendered_operator_count": rendered_count,
        "library_operator_count": len(library.operators),
        "review_version": OPERATOR_CATALOG_REVIEW_VERSION,
    }


__all__ = [
    "OPERATOR_CATALOG_REVIEW_VERSION",
    "render_generic_operator_catalog_html",
    "write_generic_operator_catalog_html",
]
