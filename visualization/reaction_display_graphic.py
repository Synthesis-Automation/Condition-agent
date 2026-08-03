"""Render a display-only minimum reaction projection as SVG."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from io import BytesIO
from typing import Any, Tuple

from reactive_taxonomy.reaction_display_projection import (
    ReactionDisplayProjection,
)

from .rendering import apply_render_preset

_DRAW_LABEL_TRANSLATION = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")

try:  # pragma: no cover - public helpers exercise availability
    from rdkit import Chem
    from rdkit.Chem import Draw, rdChemReactions
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:  # pragma: no cover
    Chem = None  # type: ignore[assignment]
    Draw = None  # type: ignore[assignment]
    rdChemReactions = None  # type: ignore[assignment]
    rdMolDraw2D = None  # type: ignore[assignment]

try:  # pragma: no cover - exercised only for annotated PNG output
    from PIL import Image, ImageDraw
except ImportError:  # pragma: no cover
    Image = None  # type: ignore[assignment]
    ImageDraw = None  # type: ignore[assignment]


@dataclass(frozen=True)
class ReactionDisplayGraphic:
    """One rendered display projection and its visible annotation."""

    image_bytes: bytes
    image_format: str
    minimum_reaction_smiles: str
    annotation: str | None
    warnings: Tuple[str, ...]


def _require_rdkit() -> None:
    if any(
        value is None
        for value in (Chem, Draw, rdChemReactions, rdMolDraw2D)
    ):
        raise RuntimeError("RDKit is required for reaction display rendering.")


def _prepare_reaction(projection: ReactionDisplayProjection) -> Any:
    labels = {
        int(value.placeholder_index): str(value.display_label)
        for value in projection.substituents
        if value.placeholder_index is not None and value.display_label
    }

    def prepare_component(component: Any) -> Any:
        molecule = Chem.MolFromSmiles(component.render_smiles)
        if molecule is None:
            raise ValueError("minimum reaction component could not be parsed")
        atom_by_placeholder = {
            int(atom.GetAtomMapNum()): int(atom.GetIdx())
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0
        }
        editable = Chem.RWMol(molecule)
        for connector in component.connectors:
            port_atoms = []
            for placeholder_index in connector.placeholder_indices:
                atom_index = atom_by_placeholder.get(int(placeholder_index))
                if atom_index is None:
                    raise ValueError(
                        "hidden connector references an unknown display port"
                    )
                port_atoms.append(atom_index)
            connector_label = str(connector.display_label).translate(
                _DRAW_LABEL_TRANSLATION
            )
            if len(port_atoms) == 2:
                if editable.GetBondBetweenAtoms(*port_atoms) is not None:
                    raise ValueError(
                        "hidden connector ports already have a chemical bond"
                    )
                editable.AddBond(
                    port_atoms[0], port_atoms[1], Chem.BondType.ZERO
                )
                connector_bond = editable.GetBondBetweenAtoms(*port_atoms)
                connector_bond.SetProp("bondNote", connector_label)
            else:
                scaffold = Chem.Atom(0)
                scaffold.SetProp("atomLabel", connector_label)
                scaffold.SetProp("_displayLabel", connector_label)
                scaffold_index = editable.AddAtom(scaffold)
                for port_atom in port_atoms:
                    editable.AddBond(
                        port_atom, scaffold_index, Chem.BondType.ZERO
                    )
        molecule = editable.GetMol()
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0:
                label = labels.get(int(atom.GetAtomMapNum()), "R")
                draw_label = label.translate(_DRAW_LABEL_TRANSLATION)
                atom.SetProp("atomLabel", draw_label)
                atom.SetProp("_displayLabel", draw_label)
            atom.SetAtomMapNum(0)
        try:
            Chem.SanitizeMol(molecule)
        except Exception:
            molecule.UpdatePropertyCache(strict=False)
        return molecule

    reaction = rdChemReactions.ChemicalReaction()
    for component in projection.reactants:
        reaction.AddReactantTemplate(prepare_component(component))
    for component in projection.products:
        reaction.AddProductTemplate(prepare_component(component))
    return reaction


def _add_svg_annotation(svg: str, annotation: str) -> str:
    """Append a short topology note inside the existing SVG canvas."""
    marker = "</svg>"
    if marker not in svg:
        return svg
    text = (
        "<text x='12' y='96%' font-family='Arial, sans-serif' "
        "font-size='15px' fill='#333333'>"
        f"{escape(annotation)}</text>"
    )
    return svg.replace(marker, f"{text}\n{marker}", 1)


def _add_png_annotation(png: bytes, annotation: str) -> bytes:
    """Add a white note strip below an RDKit PNG rendering."""
    if Image is None or ImageDraw is None:
        raise RuntimeError("Pillow is required for annotated PNG output.")
    source = Image.open(BytesIO(png)).convert("RGB")
    canvas = Image.new("RGB", (source.width, source.height + 34), "white")
    canvas.paste(source, (0, 0))
    ImageDraw.Draw(canvas).text((12, source.height + 8), annotation, fill="#333333")
    output = BytesIO()
    canvas.save(output, format="PNG")
    return output.getvalue()


def build_reaction_display_graphic(
    projection: ReactionDisplayProjection,
    *,
    size: tuple[int, int] = (1100, 280),
    image_format: str = "svg",
    render_preset: str = "current",
) -> ReactionDisplayGraphic:
    """Render a minimum reaction projection with ``R`` wildcard labels."""
    if not isinstance(projection, ReactionDisplayProjection):
        raise TypeError(
            "build_reaction_display_graphic requires a "
            "ReactionDisplayProjection"
        )
    _require_rdkit()
    width, height = size
    if width <= 0 or height <= 0:
        raise ValueError("reaction display graphic size must be positive")
    normalized_format = str(image_format).casefold()
    if normalized_format not in {"png", "svg"}:
        raise ValueError("image_format must be 'png' or 'svg'")
    reaction = _prepare_reaction(projection)
    molecules = tuple(reaction.GetReactants()) + tuple(reaction.GetProducts())
    options = rdMolDraw2D.MolDrawOptions()
    apply_render_preset(
        options,
        render_preset,
        molecules=molecules,
        context="reaction_core",
    )
    panel_count = max(
        reaction.GetNumReactantTemplates()
        + reaction.GetNumProductTemplates()
        + 1,
        3,
    )
    panel_size = (max(width // panel_count, 120), max(height, 160))
    if normalized_format == "svg":
        drawing = Draw.ReactionToImage(
            reaction,
            subImgSize=panel_size,
            useSVG=True,
            drawOptions=options,
        )
        svg = str(drawing)
        if projection.annotation:
            svg = _add_svg_annotation(svg, projection.annotation)
        image_bytes = svg.encode("utf-8")
    else:
        drawing = Draw.ReactionToImage(
            reaction,
            subImgSize=panel_size,
            returnPNG=True,
            drawOptions=options,
        )
        image_bytes = bytes(drawing)
        if projection.annotation:
            image_bytes = _add_png_annotation(
                image_bytes, projection.annotation
            )
    return ReactionDisplayGraphic(
        image_bytes=image_bytes,
        image_format=normalized_format,
        minimum_reaction_smiles=projection.minimum_reaction_smiles,
        annotation=projection.annotation,
        warnings=projection.warnings,
    )


def render_reaction_display_image_bytes(
    projection: ReactionDisplayProjection,
    *,
    size: tuple[int, int] = (1100, 280),
    image_format: str = "svg",
    render_preset: str = "current",
) -> bytes:
    """Return SVG or PNG bytes for a display-minimized reaction."""
    return build_reaction_display_graphic(
        projection,
        size=size,
        image_format=image_format,
        render_preset=render_preset,
    ).image_bytes


__all__ = [
    "ReactionDisplayGraphic",
    "build_reaction_display_graphic",
    "render_reaction_display_image_bytes",
]
