from pathlib import Path

import pytest

from reactive_taxonomy import featurize_reaction
from visualization import (
    RenderStyle,
    build_reaction_core_graphic,
    load_reaction_core_graphic_definition,
    render_molecule_image,
    render_molecule_image_bytes,
    render_reaction_core_image_bytes,
    render_reaction_image_bytes,
)

MAPPED_CLICK_REACTION = (
    "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[n:7][c:8]"
    "([N:9]([CH2:10][C:11]#[CH:12])[c:31]3[n:32][c:33]4"
    "[cH:34][cH:35][c:36]([F:37])[cH:38][c:39]4[s:40]3)"
    "[s:41][c:42]2[cH:43]1."
    "[N:13]([CH2:14][C:15](=[O:16])[NH:17][c:18]1[cH:19][cH:20]"
    "[c:21]([S:22][C:23]([F:24])([F:25])[F:26])[cH:27][cH:28]1)"
    "=[N+:29]=[N-:30]"
    ">>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[n:7][c:8]"
    "([N:9]([CH2:10][c:11]3[cH:12][n:13]"
    "([CH2:14][C:15](=[O:16])[NH:17][c:18]4[cH:19][cH:20]"
    "[c:21]([S:22][C:23]([F:24])([F:25])[F:26])[cH:27][cH:28]4)"
    "[n:29][n:30]3)[c:31]3[n:32][c:33]4[cH:34][cH:35]"
    "[c:36]([F:37])[cH:38][c:39]4[s:40]3)[s:41][c:42]2[cH:43]1"
)


def test_molecule_renderer_supports_in_memory_png_and_file_output(
    tmp_path: Path,
) -> None:
    drawing = render_molecule_image_bytes("CCO", size=(240, 160))

    assert drawing.startswith(b"\x89PNG\r\n\x1a\n")

    destination = tmp_path / "nested" / "ethanol.svg"
    result = render_molecule_image(
        "CCO",
        destination,
        size=(240, 160),
        image_format="svg",
    )

    assert result == destination
    assert "<svg" in destination.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "reaction_smiles",
    (
        "CCO>>CC=O",
        "CCO>O>CC=O",
    ),
)
def test_reaction_renderer_supports_both_reaction_smiles_forms(
    reaction_smiles: str,
) -> None:
    drawing = render_reaction_image_bytes(
        reaction_smiles,
        size=(480, 180),
    )

    assert drawing.startswith(b"\x89PNG\r\n\x1a\n")


def test_reaction_renderer_supports_vector_output() -> None:
    drawing = render_reaction_image_bytes(
        "CCO>>CC=O",
        size=(480, 180),
        image_format="svg",
    )

    assert b"<svg" in drawing[:512]


def test_reaction_core_renderer_draws_click_core_with_stable_placeholders() -> (
    None
):
    analysis = featurize_reaction(MAPPED_CLICK_REACTION)

    graphic = build_reaction_core_graphic(
        analysis,
        size=(1200, 260),
        image_format="svg",
    )
    png = render_reaction_core_image_bytes(
        analysis,
        size=(1200, 260),
        image_format="png",
    )

    assert b"<svg" in graphic.image_bytes[:512]
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert graphic.definition_id == "reaction_core_graphic.v1"
    assert graphic.schema_version == "1.0"
    assert [placeholder.label for placeholder in graphic.placeholders] == [
        "R1",
        "R2",
    ]
    assert {
        placeholder.fragment_smiles for placeholder in graphic.placeholders
    } == {
        "COc1ccc2nc(N(C)c3nc4ccc(F)cc4s3)sc2c1",
        "CC(=O)Nc1ccc(SC(F)(F)F)cc1",
    }
    assert load_reaction_core_graphic_definition()[
        "continuities_replaced_by_labels"
    ] == ["retained"]


def test_reaction_core_renderer_requires_mapped_core() -> None:
    analysis = featurize_reaction("CCO>>CC=O")

    with pytest.raises(ValueError, match="ReactionCoreProjection"):
        build_reaction_core_graphic(analysis)


@pytest.mark.parametrize(
    "reaction_smiles",
    (
        "",
        "CCO",
        ">>CCO",
        "CCO>>",
        "CCO>>CCO>>CCO",
    ),
)
def test_reaction_renderer_rejects_invalid_reaction_smiles(
    reaction_smiles: str,
) -> None:
    with pytest.raises(ValueError):
        render_reaction_image_bytes(reaction_smiles)


@pytest.mark.parametrize("size", ((0, 100), (100, -1), (100.5, 100)))
def test_render_style_rejects_invalid_sizes(size) -> None:
    with pytest.raises(ValueError):
        RenderStyle(size=size).validated_size()
