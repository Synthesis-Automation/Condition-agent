"""Regressions for the display-only reaction minimization POC."""

from reactive_taxonomy import (
    build_reaction_display_projection,
    featurize_reaction,
    reaction_render_context_from_analysis,
)


def _projection(reaction_smiles: str):
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.reaction_core is not None
    return build_reaction_display_projection(
        reaction_render_context_from_analysis(analysis)
    )


def test_suzuki_retains_aromatic_rings_and_removes_spectator_methoxy() -> None:
    projection = _projection(
        "Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1"
        ">>COc1ccc(-c2ccccc2)cc1"
    )
    assert projection.minimum_reaction_smiles == (
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    assert sum(
        value.retained_aromatic_system_count for value in projection.reactants
    ) == 2
    assert sum(
        value.removed_substituent_count for value in projection.reactants
    ) == 1
    assert sum(
        value.removed_substituent_count for value in projection.products
    ) == 1


def test_nonaromatic_hydrogenation_replaces_remote_aryl_with_r() -> None:
    projection = _projection("C=CC1=CC=CC=C1>>CCc1ccccc1")
    assert projection.minimum_reaction_smiles == "*C=C>>*CC"
    assert projection.reactants[0].retained_aromatic_system_count == 0
    assert projection.reactants[0].r_group_count == 1
    assert projection.products[0].r_group_count == 1


def test_amide_formation_retains_carboxylic_acid_carbonyl_oxygen() -> None:
    projection = _projection(
        "O=C(O)c1ccccc1.Nc1ccccc1"
        ">>O=C(Nc1ccccc1)c1ccccc1"
    )
    assert projection.minimum_reaction_smiles == (
        "*C(=O)O.*N>>*NC(*)=O"
    )
    assert "C(=*)" not in projection.minimum_reaction_smiles


def test_aliphatic_amide_formation_keeps_r_c_o_oh_handle() -> None:
    projection = _projection("CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1")
    assert projection.minimum_reaction_smiles == (
        "*C(=O)O.*N>>*NC(*)=O"
    )


def test_click_reaction_keeps_new_ring_and_uses_two_r_groups() -> None:
    projection = _projection(
        "[CH:5]#[C:6][CH2:7][CH2:8][O:9][CH3:10]."
        "[CH3:1][CH2:2][N-:3][N+:4]#[N:11]"
        ">>[CH3:1][CH2:2][n:3]1[n:4][cH:5][c:6]"
        "([CH2:7][CH2:8][O:9][CH3:10])[n:11]1"
    )
    assert projection.minimum_reaction_smiles == (
        "*C#C.*[N-][N+]#N>>*c1cnn(*)n1"
    )
    assert projection.annotation is None
    assert projection.formed_ring_sizes == (5,)
    assert sum(value.r_group_count for value in projection.reactants) == 2
    assert projection.products[0].r_group_count == 2


def test_intramolecular_projection_reports_formed_ring_size() -> None:
    projection = _projection("NCCCCBr>>C1CCCN1")
    assert projection.minimum_reaction_smiles == "*CBr.*N>>*CN*"
    assert projection.reaction_scope == "intramolecular"
    assert projection.formed_ring_sizes == (5,)
    assert projection.annotation == "Intramolecular; forms a 5-membered ring"
    assert "INTRAMOLECULAR_TETHER_ABSTRACTED_IN_DISPLAY" in (
        projection.warnings
    )


def test_projection_is_invariant_to_reactant_component_order() -> None:
    first = _projection(
        "Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1"
        ">>COc1ccc(-c2ccccc2)cc1"
    )
    second = _projection(
        "c1ccc(B(O)O)cc1.Brc1ccc(OC)cc1"
        ">>COc1ccc(-c2ccccc2)cc1"
    )
    assert first.minimum_reaction_smiles == second.minimum_reaction_smiles
