"""Regression coverage for general intra/intermolecular reaction topology."""

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.reaction_reconstruction import (
    enumerate_reconstruction_assignments,
)
from reactive_taxonomy.reaction_reconstruction_rules import (
    load_reaction_reconstruction_rules,
)


INTRAMOLECULAR_CN = "NCCc1ccccc1Br>>c1ccc2c(c1)CCN2"
MAPPED_INTRAMOLECULAR_CN = (
    "[NH2:1][CH2:2][CH2:3][c:4]1[cH:5][cH:6][cH:7][cH:8]"
    "[c:9]1[Br:10]>>"
    "[NH:1]1[CH2:2][CH2:3][c:4]2[cH:5][cH:6][cH:7][cH:8][c:9]21"
)
AZIDE_CYCLOADDITION = (
    "C#Cc1ccc(Cl)cc1."
    "[N-]=[N+]=NC1COc2c(Br)cc([N+](=O)[O-])c3cccc1c23"
    ">>O=[N+]([O-])c1cc(Br)c2c3c(cccc13)"
    "C(n1cc(-c3ccc(Cl)cc3)nn1)CO2"
)
MAPPED_DIELS_ALDER = (
    "[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]"
    ">>[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:5][CH2:6]1"
)
MAPPED_TWO_PLUS_TWO = (
    "[CH2:1]=[CH2:2].[CH2:3]=[CH2:4]"
    ">>[CH2:1]1[CH2:2][CH2:3][CH2:4]1"
)


def test_reconstruction_rules_use_general_slot_relationships() -> None:
    rules = load_reaction_reconstruction_rules()

    assert rules
    assert all(
        relationship["component_relation"] == "same_or_different"
        for rule in rules
        for relationship in rule.get("slot_relationships") or ()
    )


def test_unmapped_intramolecular_cn_cyclization_is_exactly_reconstructed() -> None:
    result = featurize_reaction(INTRAMOLECULAR_CN)

    assert result.reaction_label.concise == (
        "intramolecular (5-membered ring) Ar–Br / R–NH2 → Ar–NH–R"
    )
    assert result.reaction_label.status == "exact_reconstruction"
    assert result.transformation_class == "sp2_c_n_substitution"
    assert result.reaction_topology is not None
    assert result.reaction_topology.reaction_scope == "intramolecular"
    assert result.reaction_topology.reactant_tether_distances == (4,)
    assert result.reaction_topology.formed_ring_sizes == (5,)
    assert result.reaction_topology.ring_count_delta == 1
    assert result.reaction_topology.same_component_role_groups == ()
    assert result.reaction_topology.role_component_indices == {}
    assert result.interpretation is not None
    assert result.interpretation.selected_candidate is not None
    role_components = {
        site.component_index
        for site in result.interpretation.selected_candidate.role_assignments.values()
    }
    assert role_components == {0}


@pytest.mark.parametrize(
    ("reaction", "transformation_class", "ring_size"),
    [
        ("OCCCc1ccccc1Br>>c1ccc2c(c1)CCCO2", "sp2_c_o_substitution", 6),
        ("SCCc1ccccc1Br>>c1ccc2c(c1)CCS2", "sp2_c_s_substitution", 5),
        ("NCCCCBr>>C1CCCN1", "sp3_c_n_substitution", 5),
    ],
)
def test_shared_substitution_grammars_support_intramolecular_closure(
    reaction: str, transformation_class: str, ring_size: int
) -> None:
    result = featurize_reaction(reaction)

    assert result.reaction_label.status in {
        "exact_reconstruction",
        "family_overlay",
    }
    assert result.transformation_class == transformation_class
    assert result.reaction_topology is not None
    assert result.reaction_topology.reaction_scope == "intramolecular"
    assert result.reaction_topology.formed_ring_sizes == (ring_size,)


def test_other_shared_grammars_can_opt_into_same_component_roles() -> None:
    reductive_amination = featurize_reaction("O=CCCCN>>C1CCCN1")
    lactamization = featurize_reaction("NCCCC(=O)O>>O=C1CCCN1")

    assert reductive_amination.transformation_class == (
        "carbonyl_c_n_reductive_coupling"
    )
    assert reductive_amination.named_family == "reductive_amination"
    assert reductive_amination.reaction_topology is not None
    assert reductive_amination.reaction_topology.formed_ring_sizes == (5,)
    assert lactamization.transformation_class == "amide_formation"
    assert lactamization.reaction_topology is not None
    assert lactamization.reaction_topology.formed_ring_sizes == (5,)


def test_mapped_unknown_intramolecular_reaction_gets_generic_topology_label() -> None:
    result = featurize_reaction(
        "[Br:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2-:6]>>"
        "[CH2:2]1[CH2:3][CH2:4][CH2:5][CH2:6]1"
    )

    assert result.selected_candidate is None
    assert result.reaction_label.concise == (
        "intramolecular (5-membered ring) C–C substitution"
    )
    assert result.reaction_label.status == "generic_pattern"
    assert result.reaction_topology is not None
    assert result.reaction_topology.role_component_indices == {}
    assert result.reaction_topology.reaction_scope == "intramolecular"


@pytest.mark.parametrize(
    ("reaction", "size", "elements", "formed_bonds", "aromatic", "label"),
    (
        (
            AZIDE_CYCLOADDITION,
            5,
            ("C", "C", "N", "N", "N"),
            ("C-N", "C-N"),
            True,
            "C≡C + N=N=N → aromatic 5-membered C₂N₃ ring",
        ),
        (
            MAPPED_DIELS_ALDER,
            6,
            ("C", "C", "C", "C", "C", "C"),
            ("C-C", "C-C"),
            False,
            "C=C + C=C–C=C → 6-membered C₆ ring",
        ),
        (
            MAPPED_TWO_PLUS_TWO,
            4,
            ("C", "C", "C", "C"),
            ("C-C", "C-C"),
            False,
            "C=C + C=C → 4-membered C₄ ring",
        ),
        (
            "C=C.C=[N+]([O-])C>>C1CON(C)C1",
            5,
            ("C", "C", "O", "N", "C"),
            ("C-C", "C-O"),
            False,
            "C=C + O–N=C → 5-membered C₃NO ring",
        ),
    ),
)
def test_generic_ring_observation_and_renderer_cover_cycloadditions(
    reaction: str,
    size: int,
    elements: tuple[str, ...],
    formed_bonds: tuple[str, ...],
    aromatic: bool,
    label: str,
) -> None:
    result = featurize_reaction(reaction)

    assert result.reaction_topology is not None
    assert result.reaction_topology.formed_ring_sizes == (size,)
    assert result.reaction_topology.ring_count_delta == 1
    assert len(result.reaction_topology.ring_changes) == 1
    change = result.reaction_topology.ring_changes[0]
    assert change.ring_size == size
    assert change.element_sequence == elements
    assert change.source_component_indices == (0, 1)
    assert change.formed_bond_types == formed_bonds
    assert change.aromatic_after is aromatic
    assert result.reaction_label.concise == label
    assert result.reaction_label.status == "ring_formation"
    assert result.reaction_label is not None
    assert result.reaction_label.status == "ring_formation"
    assert result.reaction_label.source == "generic_topology"
    assert "key connectivity:" in result.reaction_label.detailed
    assert "raw edits:" in result.reaction_label.detailed


def test_ring_renderer_is_reactant_order_invariant() -> None:
    left, product = AZIDE_CYCLOADDITION.split(">>", 1)
    reversed_reactants = ".".join(reversed(left.split("."))) + ">>" + product

    forward = featurize_reaction(AZIDE_CYCLOADDITION)
    reversed_result = featurize_reaction(reversed_reactants)

    assert forward.reaction_label == reversed_result.reaction_label


def test_inter_and_intramolecular_scope_have_separate_signature_tiers() -> None:
    intramolecular = featurize_reaction(INTRAMOLECULAR_CN)
    intermolecular = featurize_reaction("Brc1ccccc1.CCN>>CCNc1ccccc1")

    assert intramolecular.reaction_signature is not None
    assert intermolecular.reaction_signature is not None
    assert intermolecular.reaction_topology is not None
    assert intermolecular.reaction_topology.reaction_scope == "intermolecular"
    assert intramolecular.reaction_signature.transformation_signature_key != (
        intermolecular.reaction_signature.transformation_signature_key
    )
    assert intramolecular.reaction_signature.bond_edit_signature_key == (
        intermolecular.reaction_signature.bond_edit_signature_key
    )


def test_mapped_and_unmapped_topology_signatures_are_identical() -> None:
    unmapped = featurize_reaction(INTRAMOLECULAR_CN)
    mapped = featurize_reaction(MAPPED_INTRAMOLECULAR_CN)

    assert unmapped.reaction_signature is not None
    assert mapped.reaction_signature is not None
    assert unmapped.reaction_signature.signature_id == (
        mapped.reaction_signature.signature_id
    )
    assert mapped.reaction_signature.topology.evidence == (
        "validated_mapping_and_exact_reconstruction"
    )


def test_topology_serializes_in_analysis_and_signature() -> None:
    payload = featurize_reaction(INTRAMOLECULAR_CN).to_dict()

    assert payload["schema_version"] == "6.0"
    assert payload["reaction_topology"]["reaction_scope"] == "intramolecular"
    assert payload["reaction_topology"]["formed_ring_sizes"] == (5,)
    assert payload["reaction_signature"]["schema_version"] == "3.2"
    assert payload["reaction_topology"]["schema_version"] == "1.2"
    assert payload["reaction_signature"]["topology"] == payload["reaction_topology"]


def test_same_and_different_relationship_constraints_are_enforced() -> None:
    intramolecular = featurize_reaction(INTRAMOLECULAR_CN)
    intermolecular = featurize_reaction("Brc1ccccc1.CCN>>CCNc1ccccc1")
    base_rule = next(
        rule
        for rule in load_reaction_reconstruction_rules()
        if rule["id"] == "sp2_c_n_substitution"
    )

    same = {
        **base_rule,
        "slot_relationships": [
            {
                "slots": ["slot_1", "slot_2"],
                "component_relation": "same",
            }
        ],
    }
    different = {
        **base_rule,
        "slot_relationships": [
            {
                "slots": ["slot_1", "slot_2"],
                "component_relation": "different",
            }
        ],
    }

    assert len(enumerate_reconstruction_assignments(intramolecular.reactants, same)) == 1
    assert enumerate_reconstruction_assignments(intramolecular.reactants, different) == ()
    assert enumerate_reconstruction_assignments(intermolecular.reactants, same) == ()
    assert len(enumerate_reconstruction_assignments(intermolecular.reactants, different)) == 1
