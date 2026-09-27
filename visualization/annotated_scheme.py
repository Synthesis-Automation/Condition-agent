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
from .scheme_labels import compact_condition_labels, compact_yield_label


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
    bottom_padding: int
    font_size: int
    line_height: int
    name_characters: int
    annotation_characters: int
    character_width_em: float
    max_annotation_lines: int
    max_name_lines: int


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
    if not isinstance(values, dict) or values.get("schema_version") != "1.1" or values.get("definition_id") != "annotated_scheme.v1":
        raise ValueError("Unsupported annotated scheme definition")
    preset = values.get("molecule_preset")
    if not isinstance(preset, str) or preset not in load_render_style_definitions()["presets"]:
        raise ValueError("Invalid scheme molecule preset")
    for key in (
        "molecule_width", "molecule_height", "arrow_width", "gap", "padding", "bottom_padding",
        "font_size", "line_height", "name_characters", "annotation_characters",
        "max_annotation_lines", "max_name_lines",
    ):
        if type(values.get(key)) is not int or not 1 <= values[key] <= 2000:
            raise ValueError(f"Invalid scheme layout parameter: {key}")
    if values["line_height"] < values["font_size"] or values["arrow_width"] < 3 * values["font_size"]:
        raise ValueError("Scheme labels require sufficient line height and arrow width")
    if (type(values.get("character_width_em")) not in {int, float}
            or not 0.5 <= values["character_width_em"] <= 1):
        raise ValueError("Invalid scheme character width estimate")
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
    """Draw a forward reaction with compact supplied labels and full metadata.

    Retrosynthetic plans are displayed in synthetic direction. No structures,
    conditions, yields, or ordering are inferred. Long labels use an ellipsis;
    complete original annotations and attribution remain in the SVG description.
    Molecule panels are minimums; larger graphs grow at the preset bond scale.
    """
    if not reactants or not products:
        raise ValueError("A scheme requires explicit reactants and products")
    style = load_annotated_scheme_style()
    mw, mh = style["molecule_width"], style["molecule_height"]
    gap, pad, aw = style["gap"], style["padding"], style["arrow_width"]
    font, line = style["font_size"], style["line_height"]

    def wrap(text: str, characters: int) -> list[str]:
        return textwrap.wrap(text, width=characters, break_long_words=True, break_on_hyphens=False) or [""]

    def bounded_lines(text: str, characters: int, cap: int) -> list[str]:
        rows = wrap(text, characters) if text else []
        if len(rows) > cap:
            rows = rows[:cap]
            rows[-1] = rows[-1][:characters - 1].rstrip(" ,;.") + "…"
        return rows

    condition_text = ", ".join(compact_condition_labels(
        (item.text for item in conditions), (item.name for item in reactants),
    ))
    # Estimate Arial label width while keeping ordinary ligand formulas intact;
    # only overlong individual tokens are split. Full text remains in metadata.
    character_width = font * style["character_width_em"]
    annotation_characters = min(style["annotation_characters"], max(1, int(aw / character_width)))
    condition_lines = bounded_lines(
        condition_text, annotation_characters, style["max_annotation_lines"],
    )
    yield_text = compact_yield_label(yield_info.text if yield_info else None)
    yield_lines = bounded_lines(yield_text or "", annotation_characters, 1)
    drawings = [ET.fromstring(render_molecule_image_bytes(
        item.smiles, size=(mw, mh), image_format="svg",
        render_preset=style["molecule_preset"], expand_canvas=True,
    )) for item in (*reactants, *products)]
    widths = [float(drawing.get("width").removesuffix("px")) for drawing in drawings]
    heights = [float(drawing.get("height").removesuffix("px")) for drawing in drawings]
    names = [bounded_lines(
        item.name, min(style["name_characters"], max(1, int(widths[index] / character_width))),
        style["max_name_lines"],
    ) for index, item in enumerate((*reactants, *products))]
    mh = max(heights)
    center_y = max(pad + mh / 2, pad + font + 24 + max(0, len(condition_lines) - 1) * line)
    molecule_top = center_y - mh / 2
    bottom = max(center_y + 38 + len(yield_lines) * line,
                 molecule_top + mh + max(map(len, names)) * line + 18)
    width = pad * 2 + sum(widths) + len(widths) * gap + aw
    height = bottom + style["bottom_padding"]
    root = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "width": str(width), "height": str(height),
        "viewBox": f"0 0 {width} {height}", "role": "img", "aria-label": title,
        "data-definition": str(style["definition_id"]),
        "data-schema-version": str(style["schema_version"]),
    })
    ET.SubElement(root, "title").text = title
    ET.SubElement(root, "desc").text = (
        f"{basis.capitalize()} transformation; synthetic direction. "
        + " + ".join(f"{m.name}: {m.smiles}" for m in reactants) + " → "
        + " + ".join(f"{m.name}: {m.smiles}" for m in products) + ". "
        + (" ".join(f"{c.basis}: {c.text}" for c in conditions) if conditions else "Conditions not supplied.")
        + (f" Yield ({yield_info.basis}): {yield_info.text}." if yield_info else " Yield not supplied.")
        + " Agent-authored display; see step evidence and limitations. Not a feasibility assessment."
    )
    ET.SubElement(root, "rect", {"width": "100%", "height": "100%", "fill": "#ffffff"})

    def label(
        x: float, y: float, text: str, *, size: int = font,
        anchor: str = "middle", role: str | None = None,
    ) -> None:
        attributes = {
            "x": str(x), "y": str(y), "text-anchor": anchor,
            "font-family": "Arial, Helvetica, sans-serif", "font-size": str(size), "fill": "#242424",
        }
        if role:
            attributes["data-role"] = role
        ET.SubElement(root, "text", attributes).text = text

    def side(items: tuple[SchemeMolecule, ...], start: float, offset: int) -> None:
        x = start
        for index, item in enumerate(items):
            drawing = drawings[offset + index]
            molecule_width = widths[offset + index]
            group = ET.SubElement(root, "g", {
                "transform": f"translate({x},{center_y - heights[offset + index] / 2})", "data-role": "molecule",
            })
            # Embed vector paths rather than nested SVG viewports; this also
            # works with SVG viewers that do not implement nested viewports.
            for child in drawing:
                for node in child.iter():
                    if node.tag.startswith("{http://www.w3.org/2000/svg}"):
                        node.tag = node.tag.split("}", 1)[1]
                group.append(child)
            for row, text in enumerate(names[offset + index]):
                label(x + molecule_width / 2, molecule_top + mh + 20 + row * line, text, role="molecule-name")
            if index < len(items) - 1:
                label(x + molecule_width + gap / 2, center_y + 8, "+", size=28)
            x += molecule_width + gap

    side(reactants, pad, 0)
    arrow_x = pad + sum(widths[:len(reactants)]) + len(reactants) * gap
    side(products, arrow_x + aw + gap, len(reactants))
    for row, text in enumerate(condition_lines):
        label(arrow_x + aw / 2, center_y - 24 - (len(condition_lines) - row - 1) * line, text, role="conditions")
    ET.SubElement(root, "path", {
        "d": f"M {arrow_x} {center_y} H {arrow_x + aw - 12}",
        "stroke": "#242424", "stroke-width": "2", "fill": "none",
        "data-role": "reaction-arrow",
    })
    ET.SubElement(root, "path", {
        "d": f"M {arrow_x + aw - 15} {center_y - 6} L {arrow_x + aw} {center_y} L {arrow_x + aw - 15} {center_y + 6} Z",
        "fill": "#242424",
    })
    for row, text in enumerate(yield_lines):
        label(arrow_x + aw / 2, center_y + 34 + row * line, text, role="yield")
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)
