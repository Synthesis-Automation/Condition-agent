"""Named reaction schemes composed from declared structures, without inference.

This is a display contract, not a reaction assessment. Annotations retain their
supplied attribution; drawing an arrow does not establish a feasible synthesis.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import textwrap
from typing import TypedDict, cast
import xml.etree.ElementTree as ET

from .rendering import load_render_style_definitions, render_molecule_image_bytes


class AnnotatedSchemeStyle(TypedDict):
    """Versioned layout values, validated before any drawing is generated."""

    schema_version: str
    definition_id: str
    molecule_preset: str
    molecule_width: int
    molecule_height: int
    arrow_width: int
    gap: int
    padding: int
    font_size: int
    line_height: int
    name_characters: int
    annotation_characters: int
    max_annotation_lines: int


@dataclass(frozen=True)
class SchemeMolecule:
    """An explicitly supplied name and structure to depict unchanged."""

    name: str
    smiles: str


@dataclass(frozen=True)
class SchemeAnnotation:
    """Display text and its stated epistemic status, not independently verified."""

    text: str
    basis: str


def load_annotated_scheme_style() -> AnnotatedSchemeStyle:
    """Load validated, versioned scheme layout parameters."""
    path = Path(__file__).with_name("definitions") / "annotated_scheme.v1.json"
    values = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(values, dict) or values.get("schema_version") != "1.0" or values.get("definition_id") != "annotated_scheme.v1":
        raise ValueError("Unsupported annotated scheme definition")
    preset = values.get("molecule_preset")
    if not isinstance(preset, str) or preset not in load_render_style_definitions()["presets"]:
        raise ValueError("Invalid scheme molecule preset")
    for key in (
        "molecule_width", "molecule_height", "arrow_width", "gap", "padding",
        "font_size", "line_height", "name_characters", "annotation_characters",
        "max_annotation_lines",
    ):
        if type(values.get(key)) is not int or not 1 <= values[key] <= 2000:
            raise ValueError(f"Invalid scheme layout parameter: {key}")
    return cast(AnnotatedSchemeStyle, values)


def render_annotated_scheme_svg(
    reactants: tuple[SchemeMolecule, ...],
    products: tuple[SchemeMolecule, ...],
    *,
    title: str,
    basis: str,
    conditions: tuple[SchemeAnnotation, ...] = (),
    yield_info: SchemeAnnotation | None = None,
) -> bytes:
    """Draw a forward reaction per step, including named molecules and labels.

Retrosynthetic plans are displayed in synthetic direction. No structures,
conditions, yields, or ordering are inferred. Long annotations have an explicit
overflow notice and remain complete in the SVG description and answer details.
"""
    if not reactants or not products:
        raise ValueError("A scheme requires explicit reactants and products")
    style = load_annotated_scheme_style()
    mw, mh = style["molecule_width"], style["molecule_height"]
    gap, pad, aw = style["gap"], style["padding"], style["arrow_width"]
    font, line = style["font_size"], style["line_height"]

    def wrap(text: str, characters: int) -> list[str]:
        return textwrap.wrap(text, width=characters, break_long_words=True, break_on_hyphens=True) or [""]

    condition_lines = [
        row for item in conditions
        for row in wrap(f"{item.basis.capitalize()}: {item.text}", style["annotation_characters"])
    ] or ["Conditions not supplied"]
    cap = style["max_annotation_lines"]
    if len(condition_lines) > cap:
        condition_lines = condition_lines[:cap - 1] + ["More conditions in step details"]
    yield_lines = wrap(
        f"Yield ({yield_info.basis}): {yield_info.text}" if yield_info else "Yield not supplied",
        style["annotation_characters"],
    )
    if len(yield_lines) > cap:
        yield_lines = yield_lines[:cap - 1] + ["More yield details in step details"]
    names = [wrap(item.name, style["name_characters"]) for item in (*reactants, *products)]
    center_y = max(140, 72 + len(condition_lines) * line)
    molecule_top = center_y - mh / 2
    bottom = max(center_y + 38 + len(yield_lines) * line,
                 molecule_top + mh + max(map(len, names)) * line + 18)
    width = pad * 2 + (len(reactants) + len(products)) * (mw + gap) + aw
    height = bottom + 62
    root = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "width": str(width), "height": str(height),
        "viewBox": f"0 0 {width} {height}", "role": "img", "aria-label": title,
        "data-definition": str(style["definition_id"]),
    })
    ET.SubElement(root, "title").text = title
    ET.SubElement(root, "desc").text = (
        f"{basis.capitalize()} transformation; synthetic direction. "
        + " + ".join(f"{m.name}: {m.smiles}" for m in reactants) + " → "
        + " + ".join(f"{m.name}: {m.smiles}" for m in products) + ". "
        + " ".join(f"{c.basis}: {c.text}" for c in conditions)
        + (f" Yield ({yield_info.basis}): {yield_info.text}." if yield_info else " Yield not supplied.")
        + " Agent-authored display; see step evidence and limitations. Not a feasibility assessment."
    )
    ET.SubElement(root, "rect", {"width": "100%", "height": "100%", "fill": "#ffffff"})

    def label(x: float, y: float, text: str, *, size: int = font, anchor: str = "middle") -> None:
        ET.SubElement(root, "text", {
            "x": str(x), "y": str(y), "text-anchor": anchor,
            "font-family": "Arial, Helvetica, sans-serif", "font-size": str(size), "fill": "#242424",
        }).text = text

    def side(items: tuple[SchemeMolecule, ...], start: float, offset: int) -> None:
        for index, item in enumerate(items):
            x = start + index * (mw + gap)
            drawing = ET.fromstring(render_molecule_image_bytes(
                item.smiles, size=(mw, mh), image_format="svg", render_preset=style["molecule_preset"],
            ))
            group = ET.SubElement(root, "g", {
                "transform": f"translate({x},{molecule_top})", "data-role": "molecule",
            })
            # Embed vector paths rather than nested SVG viewports; this also
            # works with SVG viewers that do not implement nested viewports.
            for child in drawing:
                for node in child.iter():
                    if node.tag.startswith("{http://www.w3.org/2000/svg}"):
                        node.tag = node.tag.split("}", 1)[1]
                group.append(child)
            for row, text in enumerate(names[offset + index]):
                label(x + mw / 2, molecule_top + mh + 20 + row * line, text)
            if index < len(items) - 1:
                label(x + mw + gap / 2, center_y + 8, "+", size=28)

    side(reactants, pad, 0)
    arrow_x = pad + len(reactants) * (mw + gap)
    side(products, arrow_x + aw + gap, len(reactants))
    for row, text in enumerate(condition_lines):
        label(arrow_x + aw / 2, center_y - 24 - (len(condition_lines) - row - 1) * line, text)
    ET.SubElement(root, "path", {
        "d": f"M {arrow_x} {center_y} H {arrow_x + aw - 12}",
        "stroke": "#242424", "stroke-width": "2", "fill": "none",
    })
    ET.SubElement(root, "path", {
        "d": f"M {arrow_x + aw - 15} {center_y - 6} L {arrow_x + aw} {center_y} L {arrow_x + aw - 15} {center_y + 6} Z",
        "fill": "#242424",
    })
    for row, text in enumerate(yield_lines):
        label(arrow_x + aw / 2, center_y + 34 + row * line, text)
    label(pad, height - 24, f"{basis.capitalize()} transformation · agent-authored · see step evidence and limitations", size=14, anchor="start")
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)
