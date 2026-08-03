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


def test_multisite_projection_renders_nonchemical_dashed_connector() -> None:
    projection = _projection(
        "Cc1c([N+](=O)[O-])cnc2c1c("
        "C1=CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
        ">>Cc1c(N)cnc2c1c("
        "C1CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
    )
    graphic = build_reaction_display_graphic(projection)

    assert graphic.image_bytes.count(b"stroke-dasharray") == 2


def test_true_reactant_components_remain_separate_without_connector() -> None:
    projection = _projection("CCBr.N>>CCN")
    graphic = build_reaction_display_graphic(projection)

    assert projection.connectors == ()
    assert b"stroke-dasharray" not in graphic.image_bytes


def test_dearomatization_with_exocyclic_carbonyl_renders() -> None:
    projection = _projection(
        "CCOC(=O)c1ccc("
        "C#Cc2cnc3nc(NC(=O)C(C)(C)C)[nH]c(=O)c3c2)s1"
        ">>CCOC(=O)c1ccc("
        "CCC2CNc3nc(NC(=O)C(C)(C)C)[nH]c(=O)c3C2)s1"
    )
    graphic = build_reaction_display_graphic(projection)

    assert graphic.image_bytes.startswith(b"<?xml")
    assert b"<svg" in graphic.image_bytes
