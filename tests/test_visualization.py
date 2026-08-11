import math
from dataclasses import replace
from pathlib import Path
import re

import pytest
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from reactive_taxonomy import (
    AtomMappingProviderMetadata,
    analyze_reaction_with_external_mapping,
    featurize_reaction,
    reaction_render_context_from_analysis,
    validate_external_atom_mapping,
)
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
MAPPED_REPEATED_SUZUKI = (
    "[Br:19][c:5]1[cH:6][cH:7][c:8]([Br:20])[cH:15][cH:16]1."
    "O[B:21](O)[c:4]1[cH:3][cH:2][cH:1][cH:18][cH:17]1."
    "O[B:22](O)[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1"
    ">>[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7]"
    "[c:8](-[c:9]3[cH:10][cH:11][cH:12][cH:13][cH:14]3)"
    "[cH:15][cH:16]2)[cH:17][cH:18]1"
)
MAPPED_DECARBOXYLATIVE_COUPLING = (
    "O=[C:15](O)[CH2:2][CH3:1]."
    "[cH:3]1[cH:4][c:5]([CH3:6])[n:7][c:8]2[cH:9][cH:10]"
    "[c:11]([F:12])[cH:13][c:14]12"
    ">>[CH3:1][CH2:2][c:3]1[cH:4][c:5]([CH3:6])[n:7]"
    "[c:8]2[cH:9][cH:10][c:11]([F:12])[cH:13][c:14]12"
)
MAPPED_INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION = (
    "O=[C:19](O)[C:2](=[O:1])[c:3]1[cH:4][c:5]"
    "([C:6]([F:7])([F:8])[F:9])[cH:10][cH:11][c:12]1-"
    "[c:13]1[cH:14][cH:15][cH:16][n:17][cH:18]1"
    ">>[O:1]=[C:2]1[c:3]2[cH:4][c:5]"
    "([C:6]([F:7])([F:8])[F:9])[cH:10][cH:11][c:12]2-"
    "[c:13]2[cH:14][cH:15][cH:16][n:17][c:18]21"
)
MAPPED_INTRAMOLECULAR_CN = (
    "[NH2:1][CH2:2][CH2:3][c:4]1[cH:5][cH:6][cH:7][cH:8]"
    "[c:9]1[Br:10]>>"
    "[NH:1]1[CH2:2][CH2:3][c:4]2[cH:5][cH:6][cH:7][cH:8]"
    "[c:9]21"
)
REPEATED_SUZUKI = (
    "Brc1ccc(Br)cc1.OB(O)c1ccccc1.OB(O)c1ccccc1"
    ">>c1ccc(-c2ccc(-c3ccccc3)cc2)cc1"
)
INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION = (
    "O=C(O)C(=O)c1cc(C(F)(F)F)ccc1-c1cccnc1"
    ">>O=C1c2cc(C(F)(F)F)ccc2-c2cccnc21"
)


class _RepeatedSuzukiMappingProvider:
    metadata = AtomMappingProviderMetadata(
        provider_id="rxnmapper",
        provider_version="test",
        model_id="fixture",
        model_sha256="fixture",
    )

    def map_reactions(self, reactions):
        results = []
        for reaction in reactions:
            result = validate_external_atom_mapping(
                reaction,
                MAPPED_REPEATED_SUZUKI,
                provider_metadata=self.metadata,
                mapper_confidence=0.99,
            )
            assert result.normalization is not None
            edits = tuple(
                edit
                for edit in result.normalization.edits
                if tuple(
                    sorted(
                        (
                            edit.atom_1.element,
                            edit.atom_2.element if edit.atom_2 else "H",
                        )
                    )
                )
                in {("B", "C"), ("Br", "C"), ("C", "C")}
            )
            results.append(
                replace(
                    result,
                    normalization=replace(result.normalization, edits=edits),
                )
            )
        return tuple(results)


class _IntramolecularCyclizationMappingProvider:
    metadata = AtomMappingProviderMetadata(
        provider_id="rxnmapper",
        provider_version="test",
        model_id="fixture",
        model_sha256="fixture",
    )

    def map_reactions(self, reactions):
        results = []
        for reaction in reactions:
            result = validate_external_atom_mapping(
                reaction,
                MAPPED_INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION,
                provider_metadata=self.metadata,
                mapper_confidence=0.99,
            )
            assert result.normalization is not None
            edits = tuple(
                edit
                for edit in result.normalization.edits
                if edit.atom_2 is None
                or {edit.atom_1.element, edit.atom_2.element} == {"C"}
            )
            results.append(
                replace(
                    result,
                    normalization=replace(result.normalization, edits=edits),
                )
            )
        return tuple(results)


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


def test_render_presets_are_versioned_and_include_drawing_styles() -> None:
    definition = load_render_style_definitions()

    assert definition["definition_id"] == "render_styles.v1"
    assert definition["schema_version"] == "1.0"
    assert available_render_presets() == (
        ("current", "Current"),
        ("short_bond", "Short bond"),
        ("web_consistent", "Web consistent"),
        ("acs_1996_compact", "ACS 1996"),
    )


def test_web_reaction_preset_preserves_ring_scale_across_complexity() -> None:
    simple = render_reaction_image_bytes(
        "c1ccccc1>>c1ccccc1",
        size=(980, 180),
        image_format="svg",
        render_preset="web_consistent",
    )
    coupling = render_reaction_image_bytes(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        size=(980, 180),
        image_format="svg",
        render_preset="web_consistent",
    )

    def longest_bond(drawing: bytes) -> float:
        coordinates = re.findall(
            rb"d='M ([0-9.]+),([0-9.]+) L ([0-9.]+),([0-9.]+)",
            drawing,
        )
        return max(
            math.dist(
                (float(values[0]), float(values[1])),
                (float(values[2]), float(values[3])),
            )
            for values in coordinates
        )

    assert longest_bond(simple) == pytest.approx(30.0, abs=1.1)
    assert longest_bond(coupling) == pytest.approx(
        longest_bond(simple),
        abs=1.1,
    )


def test_short_bond_preset_uses_85_percent_of_current_bond_length() -> None:
    current = render_molecule_image_bytes(
        "CCC",
        size=(480, 300),
        image_format="svg",
        render_preset="current",
    )
    short_bond = render_molecule_image_bytes(
        "CCC",
        size=(480, 300),
        image_format="svg",
        render_preset="short_bond",
    )

    assert _first_svg_bond_length(short_bond) == pytest.approx(
        _first_svg_bond_length(current) * 0.85,
        rel=0.002,
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


def test_molecule_renderer_hides_atom_mapping_numbers() -> None:
    mapped = render_molecule_image_bytes(
        "[CH3:7][OH:11]",
        size=(240, 160),
        image_format="svg",
    )
    unmapped = render_molecule_image_bytes(
        "[CH3][OH]",
        size=(240, 160),
        image_format="svg",
    )

    assert mapped == unmapped


def test_reaction_renderer_hides_atom_mapping_numbers() -> None:
    mapped = render_reaction_image_bytes(
        "[CH3:7][OH:11]>>[CH2:7]=[O:11]",
        size=(480, 180),
        image_format="svg",
    )
    unmapped = render_reaction_image_bytes(
        "[CH3][OH]>>[CH2]=[O]",
        size=(480, 180),
        image_format="svg",
    )

    assert mapped == unmapped


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
        reaction_render_context_from_analysis(analysis),
        size=(1200, 260),
        image_format="svg",
    )
    png = render_reaction_core_image_bytes(
        reaction_render_context_from_analysis(analysis),
        size=(1200, 260),
        image_format="png",
    )

    assert b"<svg" in graphic.image_bytes[:512]
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    assert graphic.definition_id == "reaction_core_graphic.v2.0"
    assert graphic.schema_version == "2.0"
    assert [placeholder.label for placeholder in graphic.placeholders] == [
        "Alk¹",
        "Alk²",
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
        "collapse_retained_multisite_scaffolds"
    ] is True
    assert load_reaction_core_graphic_definition()[
        "preserve_intramolecular_tethers"
    ] is True
    assert load_reaction_core_graphic_definition()[
        "render_nonretained_subgraphs_explicitly"
    ] is True
    assert load_reaction_core_graphic_definition()[
        "explicit_retained_subgraphs"
    ] == {
        "max_heavy_atom_count": 1,
        "remote_classes": ["heteroatom"],
    }
    assert "remote_class_labels" not in load_reaction_core_graphic_definition()


def test_reaction_core_renderer_supports_compact_acs_style() -> None:
    analysis = featurize_reaction(MAPPED_CLICK_REACTION)

    graphic = build_reaction_core_graphic(
        reaction_render_context_from_analysis(analysis),
        size=(1200, 260),
        image_format="svg",
        render_preset="acs_1996_compact",
    )

    assert b"<svg" in graphic.image_bytes[:512]
    assert b"stroke-width:0.6px" in graphic.image_bytes


def test_reaction_core_renderer_supports_multi_port_ring_boundaries() -> None:
    analysis = featurize_reaction(MAPPED_RING_ACYLATION)

    graphic = build_reaction_core_graphic(
        reaction_render_context_from_analysis(analysis),
        size=(1200, 260),
        image_format="svg",
    )
    png = render_reaction_core_image_bytes(
        reaction_render_context_from_analysis(analysis),
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


def test_reaction_core_renderer_reunites_repeated_suzuki_scaffold() -> None:
    analysis = featurize_reaction(MAPPED_REPEATED_SUZUKI)

    graphic = build_reaction_core_graphic(
        reaction_render_context_from_analysis(analysis),
        size=(1200, 260),
        image_format="svg",
    )

    assert b"<svg" in graphic.image_bytes[:512]
    assert [placeholder.label for placeholder in graphic.placeholders] == [
        "Ar¹",
        "Ar²",
        "Ar³",
    ]
    central_scaffolds = tuple(
        placeholder
        for placeholder in graphic.placeholders
        if placeholder.attachment_port_count == 2
        and placeholder.fragment_smiles == "c1ccccc1"
    )
    assert len(central_scaffolds) == 1


def test_renderer_reunites_externally_mapped_repeated_suzuki_scaffold() -> None:
    base = featurize_reaction(REPEATED_SUZUKI)
    assessment = analyze_reaction_with_external_mapping(
        REPEATED_SUZUKI,
        _RepeatedSuzukiMappingProvider(),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status == "external_mapping_only"
    graphic = build_reaction_core_graphic(
        reaction_render_context_from_analysis(assessment.analysis),
        size=(1200, 260),
        image_format="svg",
    )

    assert [placeholder.label for placeholder in graphic.placeholders] == [
        "Ar¹",
        "Ar²",
        "Ar³",
    ]
    assert sum(
        placeholder.fragment_smiles == "c1ccccc1"
        for placeholder in graphic.placeholders
    ) == 1


def test_renderer_keeps_departing_carboxyl_oxygens_explicit() -> None:
    from visualization.reaction_core_graphic import (
        _build_side_molecules,
        _placeholder_assignments,
    )

    analysis = featurize_reaction(MAPPED_DECARBOXYLATIVE_COUPLING)
    context = reaction_render_context_from_analysis(analysis)
    assignments, _, collapses = _placeholder_assignments(context)
    reactants = _build_side_molecules(
        context,
        side="reactant",
        assignments=assignments,
        scaffold_collapses=collapses,
    )

    assert any(
        sum(atom.GetAtomicNum() == 8 for atom in molecule.GetAtoms()) == 2
        for molecule in reactants
    )
    graphic = build_reaction_core_graphic(
        context,
        size=(1200, 260),
        image_format="svg",
    )
    assert b"<svg" in graphic.image_bytes[:512]


@pytest.mark.parametrize(
    "reaction_smiles",
    (
        MAPPED_INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION,
        MAPPED_INTRAMOLECULAR_CN,
    ),
)
def test_renderer_preserves_intramolecular_ring_tethers(
    reaction_smiles: str,
) -> None:
    from visualization.reaction_core_graphic import (
        _build_side_molecules,
        _placeholder_assignments,
    )

    analysis = featurize_reaction(reaction_smiles)
    context = reaction_render_context_from_analysis(analysis)
    assignments, placeholders, collapses = _placeholder_assignments(context)

    assert analysis.reaction_topology is not None
    assert analysis.reaction_topology.reaction_scope == "intramolecular"
    assert analysis.reaction_topology.ring_count_delta == 1
    assert placeholders == ()
    assert collapses == ()

    for side, components in (
        ("reactant", analysis.reactants),
        ("product", analysis.products),
    ):
        reduced = _build_side_molecules(
            context,
            side=side,
            assignments=assignments,
            scaffold_collapses=collapses,
        )
        assert len(reduced) == len(components) == 1
        expected = Chem.MolFromSmiles(components[0].input_smiles)
        assert expected is not None
        assert len(Chem.GetMolFrags(reduced[0])) == 1
        assert reduced[0].GetNumAtoms() == expected.GetNumAtoms()
        assert reduced[0].GetNumBonds() == expected.GetNumBonds()

    graphic = build_reaction_core_graphic(
        context,
        size=(1200, 260),
        image_format="svg",
    )
    assert b"<svg" in graphic.image_bytes[:512]


def test_external_mapping_preserves_intramolecular_core_topology() -> None:
    base = featurize_reaction(
        INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION
    )
    assessment = analyze_reaction_with_external_mapping(
        INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION,
        _IntramolecularCyclizationMappingProvider(),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status == "external_mapping_internal_consensus"
    core = assessment.analysis.reaction_core
    assert core is not None
    assert not hasattr(core, "abstraction")
    graphic = build_reaction_core_graphic(
        reaction_render_context_from_analysis(assessment.analysis),
        size=(1200, 260),
        image_format="svg",
    )
    assert graphic.placeholders == ()
    assert b"<svg" in graphic.image_bytes[:512]


def test_reaction_core_renderer_accepts_inferred_minimum_core() -> None:
    analysis = featurize_reaction("CCO>>CC=O")

    assert analysis.reaction_core is not None
    graphic = build_reaction_core_graphic(
        reaction_render_context_from_analysis(analysis)
    )
    assert graphic.image_bytes


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
