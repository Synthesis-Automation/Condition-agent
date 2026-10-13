"""Chemical reaction panels connected by explicitly declared route dependencies.

This renderer lays out a directed acyclic graph. It does not infer reactions,
atom correspondence, feasibility, or dependencies from molecular similarity.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import textwrap
from typing import Sequence
import xml.etree.ElementTree as ET

from .annotated_scheme import (
    SchemeAnnotation, SchemeMolecule, render_annotated_scheme_svg,
)


@dataclass(frozen=True)
class BranchedRouteStep:
    """A supplied reaction and its parents, sharing explicit molecule identifiers."""

    id: str
    reactants: tuple[SchemeMolecule, ...]
    products: tuple[SchemeMolecule, ...]
    after_step_ids: tuple[str, ...] = ()
    title: str = "Reaction step"
    basis: str = "supplied"
    conditions: tuple[SchemeAnnotation, ...] = ()
    reagents: tuple[SchemeAnnotation, ...] = ()
    yield_info: SchemeAnnotation | None = None


@dataclass(frozen=True)
class BranchedRouteStyle:
    """Versioned dimensions for panels and dependency connector lanes."""

    schema_version: str
    definition_id: str
    padding: int
    column_gap: int
    row_gap: int
    lane_gap: int
    font_size: int
    line_height: int
    title_characters: int
    arrow_head_length: int
    arrow_head_half_width: int
    stroke_width: int


def load_branched_route_style() -> BranchedRouteStyle:
    """Load validated layout definitions without executable configuration."""
    path = Path(__file__).with_name("definitions") / "branched_route_scheme.v1.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    fields = BranchedRouteStyle.__dataclass_fields__
    if not isinstance(raw, dict) or set(raw) != set(fields):
        raise ValueError("Invalid branched route definition fields")
    if (raw["schema_version"] != "1.0"
            or raw["definition_id"] != "branched_route_scheme.v1"):
        raise ValueError("Unsupported branched route definition")
    for key in fields:
        if key not in {"schema_version", "definition_id"}:
            if type(raw[key]) is not int or not 1 <= raw[key] <= 2000:
                raise ValueError(f"Invalid branched route geometry: {key}")
    if (raw["line_height"] < raw["font_size"]
            or raw["lane_gap"] < 2 * raw["arrow_head_half_width"]
            or raw["row_gap"] < 2 * raw["arrow_head_length"]):
        raise ValueError("Insufficient branched route label or connector spacing")
    return BranchedRouteStyle(**raw)


def _layers(steps: tuple[BranchedRouteStep, ...]) -> dict[str, int]:
    """Validate supplied identifiers and compute stable topological depths."""
    by_id = {step.id: step for step in steps}
    if len(by_id) != len(steps) or any(not step.id for step in steps):
        raise ValueError("Route step IDs must be nonempty and unique")
    molecules: dict[str, str] = {}
    for step in steps:
        if not step.reactants or not step.products:
            raise ValueError("Route steps require explicit reactants and products")
        for molecule in (*step.reactants, *step.products):
            if molecule.identifier:
                previous = molecules.setdefault(molecule.identifier, molecule.smiles)
                if previous != molecule.smiles:
                    raise ValueError("A molecule ID has inconsistent supplied structures")
        if len(set(step.after_step_ids)) != len(step.after_step_ids):
            raise ValueError("Duplicate route dependencies")
        for parent in step.after_step_ids:
            if parent not in by_id:
                raise ValueError(f"Unknown route dependency: {parent}")
            produced = {m.identifier for m in by_id[parent].products if m.identifier}
            consumed = {m.identifier for m in step.reactants if m.identifier}
            if not produced.intersection(consumed):
                raise ValueError("Dependencies must share an explicit intermediate ID")
    depths: dict[str, int] = {}
    while len(depths) < len(steps):
        ready = [step for step in steps if step.id not in depths
                 and all(parent in depths for parent in step.after_step_ids)]
        if not ready:
            raise ValueError("Route dependencies contain a cycle")
        for step in ready:
            depths[step.id] = max((depths[p] + 1 for p in step.after_step_ids), default=0)
    return depths


def render_branched_route_scheme_svg(
    steps: Sequence[BranchedRouteStep],
    *,
    title: str = "Branched reaction route",
    reagents_only: bool = False,
) -> bytes:
    """Draw supplied reactions in top-to-bottom layers with dependency arrows.

    Parallel steps occupy separate panels; each panel retains all its reactants
    and products. Intermediates repeat at consuming steps so convergent inputs
    remain explicit. Long dependencies use exterior lanes, never crossing a
    reaction panel. Supplied order breaks ties and defines displayed step numbers.
    Procedure annotations and yields remain in metadata when reagents_only is set.
    No input objects are mutated and no scientific assessment is performed.
    """
    supplied = tuple(steps)
    if not supplied or any(not isinstance(s, BranchedRouteStep) for s in supplied):
        raise ValueError("Supply a nonempty sequence of BranchedRouteStep objects")
    depths = _layers(supplied)
    style = load_branched_route_style()
    drawings: dict[str, ET.Element] = {}
    sizes: dict[str, tuple[float, float]] = {}
    headers: dict[str, list[str]] = {}
    for number, step in enumerate(supplied, 1):
        drawing = ET.fromstring(render_annotated_scheme_svg(
            step.reactants, step.products, title=step.title, basis=step.basis,
            conditions=step.reagents + (() if reagents_only else step.conditions),
            yield_info=None if reagents_only else step.yield_info,
        ))
        width = float(drawing.get("width"))
        lines = textwrap.wrap(
            f"{number}. {step.title} ({step.basis})", style.title_characters,
            break_long_words=True, break_on_hyphens=False,
        )
        width = max(width, max(map(len, lines)) * style.font_size + 2 * style.padding)
        height = float(drawing.get("height")) + len(lines) * style.line_height
        drawings[step.id], sizes[step.id], headers[step.id] = drawing, (width, height), lines

    rows = [[s for s in supplied if depths[s.id] == depth]
            for depth in range(max(depths.values()) + 1)]
    edges = [(parent, step.id) for step in supplied for parent in step.after_step_ids]
    long_edges = [edge for edge in edges if depths[edge[1]] > depths[edge[0]] + 1]
    row_widths = [sum(sizes[s.id][0] for s in row)
                  + (len(row) - 1) * style.column_gap for row in rows]
    content_width = max(row_widths)
    centers: dict[str, float] = {}
    for row, row_width in zip(rows, row_widths):
        x = (content_width - row_width) / 2
        for step in row:
            centers[step.id] = x + sizes[step.id][0] / 2
            x += sizes[step.id][0] + style.column_gap
    left_edges = [e for e in long_edges if centers[e[0]] <= centers[e[1]]]
    right_edges = [e for e in long_edges if e not in left_edges]
    left_gutter = (len(left_edges) + 1) * style.lane_gap if left_edges else 0
    right_gutter = (len(right_edges) + 1) * style.lane_gap if right_edges else 0
    positions: dict[str, tuple[float, float]] = {}
    # Each boundary reserves a separate horizontal lane for every endpoint.
    # This also leaves room for multiple incoming arrows at a convergence.
    boundary_edges = [[e for e in edges if depths[e[0]] == depth
                       or (e in long_edges and depths[e[1]] - 1 == depth)]
                      for depth in range(len(rows))]
    gap = max(style.row_gap, (max(map(len, boundary_edges)) + 2) * style.lane_gap)
    row_bottoms: list[float] = []
    y = float(style.padding)
    for row, row_width in zip(rows, row_widths):
        x = style.padding + left_gutter + (content_width - row_width) / 2
        for step in row:
            positions[step.id] = (x, y)
            x += sizes[step.id][0] + style.column_gap
        row_bottoms.append(y + max(sizes[s.id][1] for s in row))
        y = row_bottoms[-1] + gap
    width = 2 * style.padding + left_gutter + content_width + right_gutter
    height = y - gap + style.padding
    root = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "width": f"{width:g}",
        "height": f"{height:g}", "viewBox": f"0 0 {width:g} {height:g}",
        "role": "img", "aria-label": title, "data-layout": "layered-dependencies",
        "data-definition": style.definition_id, "data-schema-version": style.schema_version,
    })
    ET.SubElement(root, "title").text = title
    ET.SubElement(root, "desc").text = (
        "Follow dependency arrows from top to bottom; read each reaction left to right. "
        "Parallel preparations join at their consuming step. Intermediates repeat in "
        "reaction panels. Supplied transformations, not a feasibility assessment."
    )
    ET.SubElement(root, "metadata").text = json.dumps({
        "schema_version": "branched_route_scheme.v1", "style": asdict(style),
        "steps": [asdict(s) for s in supplied], "reagents_only": reagents_only,
    }, ensure_ascii=False, sort_keys=True)
    ET.SubElement(root, "rect", {"width": "100%", "height": "100%", "fill": "white"})
    for parent, child in edges:
        px, py = positions[parent]
        cx, cy = positions[child]
        pw, ph = sizes[parent]
        cw, _ = sizes[child]
        incoming = sorted((p for p, c in edges if c == child), key=centers.get)
        outgoing = sorted((c for p, c in edges if p == parent), key=centers.get)
        start_x = px + pw * (outgoing.index(child) + 1) / (len(outgoing) + 1)
        end_x = cx + cw * (incoming.index(parent) + 1) / (len(incoming) + 1)
        start_y = py + ph
        # Route from below the tallest panel in the parent's layer.
        edge = (parent, child)
        departure = depths[parent]
        bend_y = row_bottoms[departure] + (
            boundary_edges[departure].index(edge) + 1
        ) * style.lane_gap
        path = f"M {start_x:g} {start_y:g} V {bend_y:g}"
        if edge in long_edges:
            lane_x = (
                style.padding + left_edges.index(edge) * style.lane_gap
                if edge in left_edges else
                width - style.padding - right_edges.index(edge) * style.lane_gap
            )
            arrival = depths[child] - 1
            approach_y = row_bottoms[arrival] + (
                boundary_edges[arrival].index(edge) + 1
            ) * style.lane_gap
            path += f" H {lane_x:g} V {approach_y:g} H {end_x:g}"
        else:
            path += f" H {end_x:g}"
        path += f" V {cy - style.arrow_head_length:g}"
        group = ET.SubElement(root, "g", {
            "data-role": "dependency", "data-parent": parent, "data-child": child,
        })
        ET.SubElement(group, "path", {
            "d": path, "fill": "none", "stroke": "#668477",
            "stroke-width": str(style.stroke_width),
        })
        head, half = style.arrow_head_length, style.arrow_head_half_width
        ET.SubElement(group, "path", {
            "d": f"M {end_x - half:g} {cy - head:g} L {end_x:g} {cy:g} "
                 f"L {end_x + half:g} {cy - head:g} Z", "fill": "#668477",
        })
    for step in supplied:
        x, y = positions[step.id]
        w, h = sizes[step.id]
        panel = ET.SubElement(root, "g", {
            "transform": f"translate({x:g},{y:g})", "data-role": "route-step",
            "data-step-id": step.id, "data-depth": str(depths[step.id]),
            "data-width": f"{w:g}", "data-height": f"{h:g}",
        })
        ET.SubElement(panel, "rect", {
            "width": f"{w:g}", "height": f"{h:g}", "fill": "white",
            "stroke": "#a8bfae", "stroke-width": str(style.stroke_width),
        })
        for line, label in enumerate(headers[step.id]):
            ET.SubElement(panel, "text", {
                "x": str(style.padding), "y": str(style.font_size + line * style.line_height),
                "font-family": "Arial, Helvetica, sans-serif", "font-size": str(style.font_size),
                "fill": "#16342f",
            }).text = label
        drawing = drawings[step.id]
        drawing_x = (w - float(drawing.get("width"))) / 2
        drawing_y = len(headers[step.id]) * style.line_height
        content = ET.SubElement(panel, "g", {
            "transform": f"translate({drawing_x:g},{drawing_y:g})",
            "data-role": "reaction-panel",
        })
        # Inline vector paths for viewers without nested-viewport support.
        # Resolve percentage backgrounds against the original drawing canvas.
        for child in drawing:
            for node in child.iter():
                node.tag = node.tag.removeprefix("{http://www.w3.org/2000/svg}")
                for dimension in ("width", "height"):
                    if node.get(dimension) == "100%":
                        node.set(dimension, drawing.get(dimension))
            content.append(child)
        # Draw the border last so embedded white molecular canvases cannot erase it.
        ET.SubElement(panel, "rect", {
            "width": f"{w:g}", "height": f"{h:g}", "fill": "none",
            "stroke": "#a8bfae", "stroke-width": str(style.stroke_width),
        })
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)
