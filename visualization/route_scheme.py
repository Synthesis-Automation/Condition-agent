"""Three-column, row-major SVG schemes for declared linear reaction routes."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import textwrap
from typing import Sequence
import xml.etree.ElementTree as ET

from reactive_taxonomy.reaction_sequence import connect_reaction_sequence

from .annotated_scheme import SchemeAnnotation
from .rendering import load_render_style_definitions, render_molecule_image_bytes


@dataclass(frozen=True)
class RouteSchemeStep:
    """One supplied reaction; annotations are displayed without inference.

    product_index selects a zero-based dot-separated product occurrence when
    necessary. Reagent labels are separate from procedure annotations.
    Conditions, yields and evidence basis are retained in metadata.
    """

    reaction_smiles: str
    above: tuple[SchemeAnnotation, ...] = ()
    below: tuple[SchemeAnnotation, ...] = ()
    yield_info: SchemeAnnotation | None = None
    basis: str = "supplied"
    product_index: int | None = None
    reagents: tuple[SchemeAnnotation, ...] = ()


@dataclass(frozen=True)
class RouteSchemeStyle:
    """Validated versioned geometry for a route with three slots per row."""

    schema_version: str
    definition_id: str
    columns: int
    molecule_preset: str
    minimum_column_width: int
    minimum_arrow_width: int
    minimum_row_height: int
    column_gap: int
    row_gap: int
    padding: int
    block_gap: int
    font_size: int
    line_height: int
    annotation_characters: int
    character_width_em: float
    arrow_length_scale: float
    arrow_head_length: int
    arrow_head_half_height: int
    arrow_stroke_width: int


def load_route_scheme_style() -> RouteSchemeStyle:
    """Load the route style, rejecting unsupported or unsafe geometry."""
    path = Path(__file__).with_name("definitions") / "route_scheme.v1.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or set(raw) != set(
        RouteSchemeStyle.__dataclass_fields__
    ):
        raise ValueError("Invalid route scheme definition fields")
    if raw["schema_version"] != "1.2" or raw["definition_id"] != "route_scheme.v1":
        raise ValueError("Unsupported route scheme definition")
    for key in RouteSchemeStyle.__dataclass_fields__:
        if key not in {
            "schema_version",
            "definition_id",
            "molecule_preset",
            "character_width_em",
            "arrow_length_scale",
        }:
            if type(raw[key]) is not int or not 1 <= raw[key] <= 2000:
                raise ValueError(f"Invalid route scheme geometry: {key}")
    if raw["columns"] != 3 or raw["line_height"] < raw["font_size"]:
        raise ValueError(
            "Route scheme requires three columns and sufficient line height"
        )
    if (
        type(raw["character_width_em"]) not in {int, float}
        or not 1 <= raw["character_width_em"] <= 2
    ):
        raise ValueError("Invalid route scheme character width")
    presets = load_render_style_definitions()["presets"]
    preset = raw["molecule_preset"]
    if (
        not isinstance(preset, str)
        or preset not in presets
        or not presets[preset].get("target_bond_length_pixels")
    ):
        raise ValueError(
            "Route scheme requires a molecule preset with fixed bond scale"
        )
    if (
        type(raw["arrow_length_scale"]) not in {int, float}
        or not 0.25 <= raw["arrow_length_scale"] <= 1
    ):
        raise ValueError("Invalid route scheme arrow length scale")
    if (
        min(raw["minimum_column_width"], raw["minimum_arrow_width"])
        * raw["arrow_length_scale"]
        < 4 * raw["arrow_head_length"]
    ):
        raise ValueError("Route scheme arrow is too short")
    return RouteSchemeStyle(**raw)


@dataclass
class _Block:
    node: ET.Element
    width: float
    above: float
    below: float


def _group(role: str) -> ET.Element:
    return ET.Element("g", {"data-role": role})


def _place(parent: ET.Element, block: _Block, x: float, y: float) -> None:
    block.node.set("transform", f"translate({x:g},{y:g})")
    parent.append(block.node)


def _text(text: str, style: RouteSchemeStyle, role: str) -> _Block:
    lines = [
        line
        for paragraph in text.splitlines()
        for line in (
            textwrap.wrap(
                paragraph,
                style.annotation_characters,
                break_long_words=True,
                break_on_hyphens=False,
            )
            or [""]
        )
    ]
    node = _group(role)
    width = (
        max((len(line) for line in lines), default=0)
        * style.font_size
        * style.character_width_em
    )
    for i, line in enumerate(lines):
        ET.SubElement(
            node,
            "text",
            {
                "x": f"{width / 2:g}",
                "y": str(style.font_size + i * style.line_height),
                "text-anchor": "middle",
                "font-family": "Arial, Helvetica, sans-serif",
                "font-size": str(style.font_size),
                "fill": "#111111",
            },
        ).text = line
    return _Block(node, width, 0, len(lines) * style.line_height)


def _molecules(smiles: tuple[str, ...], style: RouteSchemeStyle, role: str) -> _Block:
    node = _group(role)
    drawings = []
    for value in smiles:
        root = ET.fromstring(
            render_molecule_image_bytes(
                value,
                size=(1, 1),
                image_format="svg",
                expand_canvas=True,
                render_preset=style.molecule_preset,
            )
        )
        # Inline paths retain their original scale; no raster images or rescaling.
        for child in root.iter():
            child.tag = child.tag.removeprefix("{http://www.w3.org/2000/svg}")
        width = float(root.get("width").removesuffix("px"))
        height = float(root.get("height").removesuffix("px"))
        group = _group("molecule")
        group.set("data-smiles", value)
        group.extend(list(root))
        drawings.append(_Block(group, width, height / 2, height / 2))
    height = max(block.above for block in drawings)
    x = 0.0
    for i, block in enumerate(drawings):
        _place(node, block, x, -block.above)
        x += block.width
        if i + 1 < len(drawings):
            plus = _text("+", style, "plus")
            _place(node, plus, x + style.block_gap, -style.font_size / 2)
            x += plus.width + 2 * style.block_gap
    return _Block(node, x, height, height)


def _stack(blocks: list[_Block], style: RouteSchemeStyle) -> _Block:
    node = _group("annotations")
    width = max((block.width for block in blocks), default=0)
    y = 0.0
    for block in blocks:
        _place(node, block, (width - block.width) / 2, y + block.above)
        y += block.above + block.below + style.block_gap
    return _Block(node, width, 0, max(0, y - style.block_gap))


def _arrow(
    upper: list[_Block],
    lower: list[_Block],
    style: RouteSchemeStyle,
) -> _Block:
    node = _group("arrow-block")
    top, bottom = _stack(upper, style), _stack(lower, style)
    width = max(style.minimum_arrow_width, top.width, bottom.width)
    head, half = style.arrow_head_length, style.arrow_head_half_height
    left = width * (1 - style.arrow_length_scale) / 2
    right = width - left
    ET.SubElement(
        node,
        "path",
        {
            "d": f"M {left:g} 0 H {right - head:g}",
            "stroke": "#111111",
            "stroke-width": str(style.arrow_stroke_width),
            "fill": "none",
            "data-role": "reaction-arrow",
        },
    )
    ET.SubElement(
        node,
        "path",
        {
            "d": f"M {right - head:g} {-half} L {right:g} 0 L {right - head:g} {half} Z",
            "fill": "#111111",
        },
    )
    _place(node, top, (width - top.width) / 2, -style.block_gap - top.below)
    _place(node, bottom, (width - bottom.width) / 2, style.block_gap)
    return _Block(
        node,
        width,
        max(half, top.below + style.block_gap),
        max(half, bottom.below + style.block_gap),
    )


def render_route_scheme_svg(
    steps: Sequence[str | RouteSchemeStep],
    *,
    title: str = "Reaction route",
    main_reactant_index: int | None = None,
    reagents_only: bool = False,
) -> bytes:
    """Draw a linear route in three columns, reading left-to-right then down.

    Each intermediate is drawn once. All supplied components are retained:
    partners and agents above arrows, unselected products below. The first
    reactants share one block unless a main reactant is explicitly selected.
    This is a depiction of supplied chemistry, not a feasibility assessment.
    With reagents_only, display explicit reagent annotations and structures;
    procedure annotations and yield remain in metadata rather than on arrows.
    """
    if isinstance(steps, str) or not steps:
        raise ValueError("Supply a nonempty ordered sequence of route steps")
    supplied = tuple(
        RouteSchemeStep(step) if isinstance(step, str) else step for step in steps
    )
    if any(not isinstance(step, RouteSchemeStep) for step in supplied):
        raise ValueError("Each route step must be reaction SMILES or RouteSchemeStep")
    connected = connect_reaction_sequence(
        tuple(step.reaction_smiles for step in supplied),
        main_reactant_index=main_reactant_index,
        product_indices=tuple(step.product_index for step in supplied),
    )
    style = load_route_scheme_style()
    first = connected[0]
    start = (
        first.reactants
        if main_reactant_index is None
        else (first.reactants[main_reactant_index],)
    )
    blocks = [_molecules(start, style, "route-molecule")]
    for number, (step, connection) in enumerate(zip(supplied, connected), 1):
        partners = tuple(
            value
            for i, value in enumerate(connection.reactants)
            if i != connection.carried_reactant_index
        )
        if number == 1 and main_reactant_index is None:
            partners = ()
        upper = []
        if partners:
            upper.append(_molecules(partners, style, "partners"))
        if connection.agents:
            upper.append(_molecules(connection.agents, style, "agents"))
        upper.extend(
            _text(item.text, style, "reagents")
            for item in step.reagents
            if item.text.strip()
        )
        upper.extend(
            _text(item.text, style, "conditions-above")
            for item in step.above
            if not reagents_only and item.text.strip()
        )
        lower = [
            _text(item.text, style, "conditions-below")
            for item in step.below
            if not reagents_only and item.text.strip()
        ]
        if not reagents_only and step.yield_info and step.yield_info.text.strip():
            lower.append(_text(step.yield_info.text, style, "yield"))
        other = tuple(
            value
            for i, value in enumerate(connection.products)
            if i != connection.product_index
        )
        if other:
            lower.extend(
                (
                    _text("Other products", style, "other-products-label"),
                    _molecules(other, style, "other-products"),
                )
            )
        arrow = _arrow(upper, lower, style)
        arrow.node.set("data-step", str(number))
        blocks.extend(
            (
                arrow,
                _molecules(
                    (connection.products[connection.product_index],),
                    style,
                    "route-molecule",
                ),
            )
        )
    rows = [blocks[i : i + style.columns] for i in range(0, len(blocks), style.columns)]
    # A wrapped route puts arrows in different columns on successive rows.
    # Size each slot by its content, so a large molecule cannot widen an arrow.
    row_widths = [
        [
            block.width
            if block.node.get("data-role") == "arrow-block"
            else max(style.minimum_column_width, block.width)
            for block in row
        ]
        for row in rows
    ]
    column_width = max(value for row in row_widths for value in row)
    extents = [
        (
            max(style.minimum_row_height / 2, *(b.above for b in row)),
            max(style.minimum_row_height / 2, *(b.below for b in row)),
        )
        for row in rows
    ]
    width = (
        2 * style.padding
        + max(
            sum(widths) + (len(widths) - 1) * style.column_gap
            for widths in row_widths
        )
    )
    height = (
        2 * style.padding
        + sum(a + b for a, b in extents)
        + (len(rows) - 1) * style.row_gap
    )
    root = ET.Element(
        "svg",
        {
            "xmlns": "http://www.w3.org/2000/svg",
            "width": f"{width:g}",
            "height": f"{height:g}",
            "viewBox": f"0 0 {width:g} {height:g}",
            "role": "img",
            "aria-label": title,
            "data-definition": style.definition_id,
            "data-schema-version": style.schema_version,
            "data-layout": "three-column-row-major",
            "data-column-width": f"{column_width:g}",
            "data-spacing": "content-sized",
        },
    )
    ET.SubElement(root, "title").text = title
    description = (
        "Read left to right, then continue at the left of the next row. "
        "Supplied transformations; not a feasibility assessment. "
        "No conditions or yields are inferred. See metadata for attribution."
    )
    if reagents_only:
        description += (
            " Only reagents label arrows; full procedures and yields are in metadata."
        )
    ET.SubElement(root, "desc").text = description
    ET.SubElement(root, "metadata").text = json.dumps(
        {
            "schema_version": "route_scheme.v1",
            "style": asdict(style),
            "steps": [asdict(step) for step in supplied],
            "connections": [asdict(step) for step in connected],
            "main_reactant_index": main_reactant_index,
            "reagents_only": reagents_only,
        },
        ensure_ascii=False,
        sort_keys=True,
    )
    ET.SubElement(root, "rect", {"width": "100%", "height": "100%", "fill": "white"})
    y = float(style.padding)
    for row_index, (row, (above, below)) in enumerate(zip(rows, extents)):
        x = float(style.padding)
        for column, block in enumerate(row):
            block.node.set("data-row", str(row_index))
            block.node.set("data-column", str(column))
            block.node.set("data-width", f"{block.width:g}")
            block.node.set("data-above", f"{block.above:g}")
            block.node.set("data-below", f"{block.below:g}")
            slot_width = row_widths[row_index][column]
            _place(root, block, x + (slot_width - block.width) / 2, y + above)
            x += slot_width + style.column_gap
        y += above + below + style.row_gap
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)
