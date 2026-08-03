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
    reaction = rdChemReactions.ReactionFromSmarts(
        projection.minimum_reaction_smiles,
        useSmiles=True,
    )
    if reaction is None:
        raise ValueError("minimum reaction SMILES could not be parsed")
    molecules = tuple(reaction.GetReactants()) + tuple(reaction.GetProducts())
    for molecule in molecules:
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
            if atom.GetAtomicNum() == 0:
                atom.SetProp("atomLabel", "R")
                atom.SetProp("_displayLabel", "R")
        try:
            Chem.SanitizeMol(molecule)
        except Exception:
            molecule.UpdatePropertyCache(strict=False)
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
