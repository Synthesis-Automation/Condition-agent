"""Regressions for the display-only reaction minimization POC."""

import pytest

from reactive_taxonomy import (
    build_reaction_display_projection,
    featurize_reaction,
    reaction_render_context_from_analysis,
)
from reactive_taxonomy.reaction_display_projection import (
    load_reaction_interface_block_definitions,
)


MAPPED_INTERCOMPONENT_ANNULATION = (
    "CS(=O)(=O)O[CH2:1][CH2:2][c:3]1[c:9](Cl)"
    "[cH:8][cH:7][cH:6][cH:4]1.[NH2:28][CH3:26]>>"
    "[CH3:26][N:28]1[CH2:1][CH2:2][c:3]2[c:9]1"
    "[cH:8][cH:7][cH:6][cH:4]2"
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
    hidden = [
        value
        for value in projection.substituents
        if value.boundary_action == "aromatic_hydrogen_cap"
    ]
    assert len(hidden) == 2
    assert {value.fragment_smiles for value in hidden} == {"CO"}
    assert {
        (
            relation.positional_relation,
            relation.aromatic_ring_distance,
            relation.same_ring,
        )
        for value in hidden
        for relation in value.aromatic_relations
    } == {("para", 3, True)}


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


def test_ester_hydrolysis_retains_acyl_pi_system_from_active_oxygen() -> None:
    projection = _projection(
        "[CH3:1][O:2][C:3]([CH3:4])=[O:5]"
        ">>[OH:2][C:3]([CH3:4])=[O:5]"
    )

    assert projection.minimum_reaction_smiles == (
        "*C(=O)OC>>*C(=O)O"
    )
    assert projection.render_reaction_smiles == (
        "COC(=O)[*:1]>>O=C(O)[*:1]"
    )


def test_ordinary_ether_cleavage_does_not_gain_pi_system_context() -> None:
    projection = _projection(
        "[CH3:1][O:2][CH3:3]>>[OH:2][CH3:3]"
    )

    assert projection.minimum_reaction_smiles == "*OC>>*O"
    assert projection.render_reaction_smiles == "CO[*:1]>>O[*:1]"


def test_reaction_interface_block_registry_is_validated() -> None:
    definitions = load_reaction_interface_block_definitions()

    assert {definition["id"] for definition in definitions} == {
        "acyl_heteroatom",
        "azide",
        "isocyanate",
        "isothiocyanate",
        "nitrile",
        "phosphoryl_heteroatoms",
        "sulfonyl_heteroatom",
        "thiocyanate",
    }
    assert all(
        definition["patterns"] and definition["retained_atom_maps"]
        for definition in definitions
    )


@pytest.mark.parametrize(
    ("reaction_smiles", "minimum_reaction_smiles"),
    (
        (
            "CC(Br)N=[N+]=[N-].N>>CC(N)N=[N+]=[N-]",
            "*C(Br)N=[N+]=[N-].N>>*C(N)N=[N+]=[N-]",
        ),
        (
            "CC(Br)[N-][N+]#N.N>>CC(N)[N-][N+]#N",
            "*C(Br)[N-][N+]#N.N>>*C(N)[N-][N+]#N",
        ),
        (
            "CC(Br)SC#N.N>>CC(N)SC#N",
            "*C(Br)SC#N.N>>*C(N)SC#N",
        ),
        (
            "CC(Br)N=C=S.N>>CC(N)N=C=S",
            "*C(Br)N=C=S.N>>*C(N)N=C=S",
        ),
        (
            "CC(Br)N=C=O.N>>CC(N)N=C=O",
            "*C(Br)N=C=O.N>>*C(N)N=C=O",
        ),
        (
            "CC(Br)C#N.N>>CC(N)C#N",
            "*C(Br)C#N.N>>*C(N)C#N",
        ),
        (
            "CC(Br)S(=O)(=O)Cl.N>>CC(N)S(=O)(=O)Cl",
            "*C(Br)S(=O)(=O)Cl.N>>*C(N)S(=O)(=O)Cl",
        ),
        (
            "CC(Br)P(=O)(O)Cl.N>>CC(N)P(=O)(O)Cl",
            "*C(Br)P(=O)(O)Cl.N>>*C(N)P(=O)(O)Cl",
        ),
    ),
)
def test_complete_interface_blocks_are_retained_at_core_boundary(
    reaction_smiles: str,
    minimum_reaction_smiles: str,
) -> None:
    projection = _projection(reaction_smiles)

    assert projection.minimum_reaction_smiles == minimum_reaction_smiles


def test_remote_nitrile_does_not_expand_through_hidden_carbon_chain() -> None:
    projection = _projection("BrCCC#N.N>>NCCC#N")

    assert projection.minimum_reaction_smiles == "*CBr.N>>*CN"


def test_thiocyanate_and_isothiocyanate_remain_topologically_distinct() -> None:
    thiocyanate = _projection("CC(Br)SC#N.N>>CC(N)SC#N")
    isothiocyanate = _projection("CC(Br)N=C=S.N>>CC(N)N=C=S")

    assert "SC#N" in thiocyanate.minimum_reaction_smiles
    assert "N=C=S" in isothiocyanate.minimum_reaction_smiles
    assert (
        thiocyanate.minimum_reaction_smiles
        != isothiocyanate.minimum_reaction_smiles
    )


def test_multisite_amide_formation_and_ester_hydrolysis_keep_both_acyls() -> None:
    projection = _projection(
        "COC(=O)[C@H](CSCCCC(=O)O)"
        "NC(=O)OCC1c2ccccc2-c2ccccc21.Nc1ccccc1"
        ">>O=C(CCCSC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc21)"
        "C(=O)O)Nc1ccccc1"
    )

    assert projection.minimum_reaction_smiles == (
        "*C(=O)O.*C(=O)OC.*N>>*C(=O)O.*NC(*)=O"
    )
    assert projection.render_reaction_smiles == (
        "COC(=O)[*:1].O=C(O)[*:2].N[*:3]"
        ">>O=C(N[*:3])[*:1].O=C(O)[*:2]"
    )
    assert {
        (connector.side, connector.port_display_labels)
        for connector in projection.connectors
    } == {
        ("reactant", ("R¹", "R²")),
        ("product", ("R¹", "R²")),
    }


@pytest.mark.parametrize(
    ("reaction_smiles", "minimum_reaction_smiles"),
    (
        ("CCBr.N>>CCN", "*CBr.N>>*CN"),
        ("CC=CBr.N>>CC=CN", "*C=CBr.N>>*C=CN"),
        ("CC#CBr.N>>CC#CN", "*C#CBr.N>>*C#CN"),
        ("CC(=O)Cl.CN>>CC(=O)NC", "*C(=O)Cl.*N>>*NC(*)=O"),
    ),
)
def test_reaction_interface_closure_follows_local_bonding(
    reaction_smiles: str,
    minimum_reaction_smiles: str,
) -> None:
    projection = _projection(reaction_smiles)
    assert projection.minimum_reaction_smiles == minimum_reaction_smiles


def test_noncarbon_x_h_center_records_distinct_substituent_groups() -> None:
    projection = _projection(
        "Brc1ccccc1.CNC>>CN(C)c1ccccc1"
    )
    assert projection.minimum_reaction_smiles == (
        "*N*.Brc1ccccc1>>*N(*)c1ccccc1"
    )
    assert projection.render_reaction_smiles == (
        "N([*:1])[*:2].Brc1ccccc1"
        ">>c1ccc(N([*:1])[*:2])cc1"
    )
    nitrogen_substituents = [
        value
        for value in projection.substituents
        if value.center_element == "N" and value.display_label
    ]
    assert {value.display_label for value in nitrogen_substituents} == {
        "R¹",
        "R²",
    }
    assert all(value.fragment_smiles == "C" for value in nitrogen_substituents)


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


def test_intercomponent_annulation_retains_formed_ring_path() -> None:
    projection = _projection(MAPPED_INTERCOMPONENT_ANNULATION)

    assert projection.definition_id == "reaction_display_projection.v1.9"
    assert projection.reaction_scope == "intermolecular"
    assert projection.formed_ring_sizes == (5,)
    assert projection.minimum_reaction_smiles == (
        "*COS(C)(=O)=O.*c1ccccc1Cl.*N>>*CN(*)c1ccccc1*"
    )
    assert {
        (connector.side, connector.port_display_labels)
        for connector in projection.connectors
    } == {
        ("reactant", ("R¹", "R²")),
        ("product", ("R¹", "R²")),
    }
    ring_path_ports = {
        (value.side, value.center_atom_map_number): (
            value.display_label,
            value.remote_class,
        )
        for value in projection.substituents
        if value.center_atom_map_number in {1, 3}
        and value.attachment_atom_map_number == 2
    }
    assert ring_path_ports == {
        ("reactant", 1): ("R¹", "alkyl"),
        ("reactant", 3): ("R²", "alkyl"),
        ("product", 1): ("R¹", "ring_aliphatic"),
        ("product", 3): ("R²", "ring_aliphatic"),
    }


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
    assert first.render_reaction_smiles == second.render_reaction_smiles


def test_indexed_r_groups_are_invariant_to_reactant_order() -> None:
    first = _projection(
        "O=C(O)c1ccccc1.Nc1ccccc1"
        ">>O=C(Nc1ccccc1)c1ccccc1"
    )
    second = _projection(
        "Nc1ccccc1.O=C(O)c1ccccc1"
        ">>O=C(Nc1ccccc1)c1ccccc1"
    )
    assert first.render_reaction_smiles == second.render_reaction_smiles
    assert first.render_reaction_smiles == (
        "O=C(O)[*:1].N[*:2]>>O=C(N[*:2])[*:1]"
    )


def test_multisite_projection_records_shortest_hidden_connector() -> None:
    projection = _projection(
        "Cc1c([N+](=O)[O-])cnc2c1c("
        "C1=CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
        ">>Cc1c(N)cnc2c1c("
        "C1CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
    )

    assert len(projection.reactants) == 1
    assert len(projection.products) == 1
    assert len(projection.connectors) == 2
    assert {
        (value.side, value.connector_id, value.port_display_labels)
        for value in projection.connectors
    } == {
        ("reactant", "S1", ("R³", "R⁴")),
        ("product", "S1", ("R³", "R⁴")),
    }
    assert all(
        value.shortest_path_bond_count == 3
        and len(value.shortest_path_atom_indices) == 4
        for value in projection.connectors
    )
    assert "HIDDEN_CONNECTOR_ABSTRACTED_IN_DISPLAY" in projection.warnings


def test_terminal_omissions_do_not_create_hidden_connectors() -> None:
    projection = _projection("C=CC1=CC=CC=C1>>CCc1ccccc1")
    assert projection.connectors == ()


def test_r_label_uses_correspondence_when_stereo_smiles_text_changes() -> None:
    scaffold = (
        "CC1=C(c2ccc(OC3CCCCO3)cc2)C(c2ccc(I)cc2)"
        "Oc2ccc(OC3CCCCO3)cc21"
    )
    alcohol = "C[C@@H]1CCN([C@@H](C)CO)C1"
    product = (
        "CC1=C(c2ccc(OC3CCCCO3)cc2)"
        "C(c2ccc(OC[C@H](C)N3CC[C@@H](C)C3)cc2)"
        "Oc2ccc(OC3CCCCO3)cc21"
    )
    first = _projection(f"{scaffold}.{alcohol}>>{product}")
    reordered = _projection(f"{alcohol}.{scaffold}>>{product}")

    assert first.render_reaction_smiles == (
        "O[*:1].Ic1ccccc1>>c1ccc(O[*:1])cc1"
    )
    assert reordered.render_reaction_smiles == first.render_reaction_smiles
    labeled = [
        value for value in first.substituents if value.display_label is not None
    ]
    assert {value.display_label for value in labeled} == {"R"}
    assert len({value.atom_map_numbers for value in labeled}) == 1
    assert labeled[0].atom_map_numbers


def test_aromatic_valence_completion_retains_exocyclic_carbonyl() -> None:
    projection = _projection(
        "CCOC(=O)c1ccc("
        "C#Cc2cnc3nc(NC(=O)C(C)(C)C)[nH]c(=O)c3c2)s1"
        ">>CCOC(=O)c1ccc("
        "CCC2CNc3nc(NC(=O)C(C)(C)C)[nH]c(=O)c3C2)s1"
    )

    assert projection.minimum_reaction_smiles == (
        "*C#Cc1cnc2nc[nH]c(=O)c2c1"
        ">>*CCC1CNc2nc[nH]c(=O)c2C1"
    )
    assert "c(=O)" in projection.reactants[0].display_smiles
    assert "c(=O)" in projection.products[0].display_smiles


def test_fragment_projection_ignores_stereobonds_outside_omitted_fragment() -> None:
    """Parent stereobonds outside atomsToUse must not break serialization."""
    projection = _projection(
        "CC(C)(C)O[C:18](=O)[N:17]1[C@@H:8]("
        "[CH2:7][CH2:6][C:5](=[O:19])[CH2:4]/[CH:3]=[CH:2]/[CH3:1])"
        "[CH2:9][CH2:10][CH2:11][C@H:12]1[C:13](=[O:14])"
        "[CH2:15][CH3:16]>>[CH3:1]/[CH:2]=[CH:3]\\[CH:4]="
        "[C:5]1/[CH2:6][CH2:7][C@H:8]2[CH2:9][CH2:10][CH2:11]"
        "[C@@H:12]([C:13](=[O:14])[CH2:15][CH3:16])[N:17]12"
    )

    assert projection.minimum_reaction_smiles == (
        "*CC(*)=O.*N(*)C(=O)OC(C)(C)C>>*C=C(*)N(*)*"
    )


def test_hofmann_rearrangement_has_no_false_three_membered_ring_note() -> None:
    projection = _projection(
        "NC(=O)CC1(CC(=O)O)CCCCC1>>NCC1(CC(=O)O)CCCCC1"
    )

    assert projection.minimum_reaction_smiles == "*CC(N)=O>>*CN"
    assert projection.formed_ring_sizes == ()
    assert projection.annotation is None
    assert "INTRAMOLECULAR_TETHER_ABSTRACTED_IN_DISPLAY" not in (
        projection.warnings
    )
