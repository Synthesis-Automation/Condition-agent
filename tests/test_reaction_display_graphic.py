"""SVG rendering regressions for display-minimized reactions."""

from reactive_taxonomy import (
    build_reaction_display_projection,
    featurize_reaction,
    reaction_render_context_from_analysis,
)
from visualization import build_reaction_display_graphic


def _projection(reaction_smiles: str):
    analysis = featurize_reaction(reaction_smiles)
    return build_reaction_display_projection(
        reaction_render_context_from_analysis(analysis)
    )


def test_display_projection_renders_as_svg() -> None:
    projection = _projection("C=CC1=CC=CC=C1>>CCc1ccccc1")
    graphic = build_reaction_display_graphic(projection)
    assert graphic.image_format == "svg"
    assert graphic.image_bytes.startswith(b"<?xml")
    assert b"<svg" in graphic.image_bytes
    assert graphic.minimum_reaction_smiles == "*C=C>>*CC"


def test_display_projection_renders_as_png() -> None:
    projection = _projection("C=CC1=CC=CC=C1>>CCc1ccccc1")
    graphic = build_reaction_display_graphic(
        projection,
        image_format="png",
    )
    assert graphic.image_format == "png"
    assert graphic.image_bytes.startswith(b"\x89PNG\r\n\x1a\n")


def test_intramolecular_note_is_visible_in_svg() -> None:
    projection = _projection("NCCCCBr>>C1CCCN1")
    graphic = build_reaction_display_graphic(projection)
    assert graphic.annotation == "Intramolecular; forms a 5-membered ring"
    assert b"Intramolecular; forms a 5-membered ring" in graphic.image_bytes
