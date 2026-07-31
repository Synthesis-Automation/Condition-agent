import math
from pathlib import Path
import re

import pytest
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from reactive_taxonomy import featurize_reaction
from visualization import (
    RenderStyle,
    apply_render_preset,
    available_render_presets,
    build_reaction_core_graphic,
    load_reaction_core_graphic_definition,
    load_render_style_definitions,
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
MAPPED_RING_ACYLATION = (
    "[CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[CH:8]=[O:9]."
    "[cH:10]1[n:11][cH:12][c:13](-[c:14]2[cH:15][cH:16][cH:17]"
    "[cH:18][cH:19]2)[c:20]2[cH:21][cH:22][cH:23][cH:24][c:25]12"
    ">>[CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[C:8](=[O:9])"
    "[c:10]1[n:11][cH:12][c:13](-[c:14]2[cH:15][cH:16][cH:17]"
    "[cH:18][cH:19]2)[c:20]2[cH:21][cH:22][cH:23][cH:24][c:25]12"
)


def _first_svg_bond_length(drawing: bytes) -> float:
    text = drawing.decode("utf-8")
    line = next(line for line in text.splitlines() if "class='bond-0" in line)
    match = re.search(
        r"M ([0-9.]+),([0-9.]+) L ([0-9.]+),([0-9.]+)",
        line,
    )
    assert match is not None
    coordinates = tuple(float(value) for value in match.groups())
    return math.dist(coordinates[:2], coordinates[2:])


def test_render_presets_are_versioned_and_include_compact_acs() -> None:
    definition = load_render_style_definitions()

    assert definition["definition_id"] == "render_styles.v1"
    assert definition["schema_version"] == "1.0"
    assert available_render_presets() == (
        ("current", "Current"),
        ("acs_1996_compact", "ACS 1996"),
    )


def test_acs_preset_uses_native_bond_length_and_line_width() -> None:
    current = render_molecule_image_bytes(
        "CCC",
        size=(480, 300),
        image_format="svg",
        render_preset="current",
    )
    compact = render_molecule_image_bytes(
        "CCC",
        size=(480, 300),
        image_format="svg",
        render_preset="acs_1996_compact",
    )

    assert _first_svg_bond_length(current) > 200.0
    assert _first_svg_bond_length(compact) == pytest.approx(14.4, abs=0.2)
    assert b"stroke-width:0.6px" in compact


def test_acs_preset_preserves_native_atom_font_to_bond_proportion() -> None:
    molecule = rdMolDraw2D.PrepareMolForDrawing(Chem.MolFromSmiles("CCl"))
    mean_bond_length = rdMolDraw2D.MeanBondLength(molecule)
    options = rdMolDraw2D.MolDrawOptions()

    apply_render_preset(
        options,
        "acs_1996_compact",
        molecules=(molecule,),
    )

    rendered_bond_length = options.fixedBondLength * mean_bond_length
    assert rendered_bond_length == pytest.approx(14.4)
    assert options.fixedFontSize == 10
    assert options.fixedFontSize / rendered_bond_length > 0.69
    assert options.multipleBondOffset == pytest.approx(0.18)


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


def test_reaction_renderer_supports_compact_acs_style() -> None:
    drawing = render_reaction_image_bytes(
        "CCO>>CC=O",
        size=(480, 180),
        image_format="svg",
        render_preset="acs_1996_compact",
    )

    assert b"<svg" in drawing[:512]
    assert b"stroke-width:0.6px" in drawing


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
    assert load_reaction_core_graphic_definition()[
        "explicit_retained_subgraphs"
    ] == {
        "max_heavy_atom_count": 1,
        "remote_classes": ["heteroatom"],
    }


def test_reaction_core_renderer_supports_compact_acs_style() -> None:
    analysis = featurize_reaction(MAPPED_CLICK_REACTION)

    graphic = build_reaction_core_graphic(
        analysis,
        size=(1200, 260),
        image_format="svg",
        render_preset="acs_1996_compact",
    )

    assert b"<svg" in graphic.image_bytes[:512]
    assert b"stroke-width:0.6px" in graphic.image_bytes


def test_reaction_core_renderer_supports_multi_port_ring_boundaries() -> None:
    analysis = featurize_reaction(MAPPED_RING_ACYLATION)

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
    placeholders = {
        placeholder.label: placeholder for placeholder in graphic.placeholders
    }
    assert set(placeholders) == {"Ar", "HetAr"}
    assert placeholders["Ar"].attachment_port_count == 1
    assert placeholders["HetAr"].attachment_port_count == 2
    assert placeholders["HetAr"].fragment_smiles == (
        "ncc(-c1ccccc1)c1ccccc1"
    )
    assert all(
        placeholder.remote_class != "heteroatom"
        for placeholder in graphic.placeholders
    )


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


def test_render_style_rejects_unknown_preset() -> None:
    with pytest.raises(ValueError, match="render_preset"):
        RenderStyle(size=(100, 100), render_preset="unknown").validated_preset()
