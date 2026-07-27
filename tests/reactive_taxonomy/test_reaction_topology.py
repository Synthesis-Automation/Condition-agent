"""Regression coverage for general intra/intermolecular reaction topology."""

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.reaction_candidates import enumerate_role_assignments
from reactive_taxonomy.reaction_grammars import load_reaction_grammars


INTRAMOLECULAR_CN = "NCCc1ccccc1Br>>c1ccc2c(c1)CCN2"
MAPPED_INTRAMOLECULAR_CN = (
    "[NH2:1][CH2:2][CH2:3][c:4]1[cH:5][cH:6][cH:7][cH:8]"
    "[c:9]1[Br:10]>>"
    "[NH:1]1[CH2:2][CH2:3][c:4]2[cH:5][cH:6][cH:7][cH:8][c:9]21"
)


def test_grammars_use_general_role_relationships() -> None:
    grammars = load_reaction_grammars()

    assert grammars
    assert all("distinct_components" not in grammar for grammar in grammars)
    assert all(
        relationship["component_relation"] == "same_or_different"
        for grammar in grammars
        for relationship in grammar.get("role_relationships") or ()
    )


def test_unmapped_intramolecular_cn_cyclization_is_exactly_reconstructed() -> None:
    result = featurize_reaction(INTRAMOLECULAR_CN)

    assert result.reaction_label == (
        "intramolecular (5-membered ring) Ar–Br / R–NH2 → Ar–NH–R"
    )
    assert result.reaction_label_status == "exact_product"
    assert result.transformation_class == "sp2_c_n_substitution"
    assert result.reaction_topology is not None
    assert result.reaction_topology.reaction_scope == "intramolecular"
    assert result.reaction_topology.reactant_tether_distances == (4,)
    assert result.reaction_topology.formed_ring_sizes == (5,)
    assert result.reaction_topology.ring_count_delta == 1
    assert result.reaction_topology.same_component_role_groups == (
        ("electrophile", "nucleophile"),
    )


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

    assert result.reaction_label_status == "exact_product"
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
    assert result.reaction_label == (
        "intramolecular (5-membered ring) C–C substitution"
    )
    assert result.reaction_label_status == "mapped_generic_pattern"
    assert result.reaction_topology is not None
    assert result.reaction_topology.role_component_indices == {}
    assert result.reaction_topology.reaction_scope == "intramolecular"


def test_inter_and_intramolecular_scope_have_separate_signature_tiers() -> None:
    intramolecular = featurize_reaction(INTRAMOLECULAR_CN)
    intermolecular = featurize_reaction(
        "Brc1ccccc1.CCN>>CCNc1ccccc1"
    )

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

    assert payload["schema_version"] == "1.8"
    assert payload["reaction_topology"]["reaction_scope"] == "intramolecular"
    assert payload["reaction_topology"]["formed_ring_sizes"] == (5,)
    assert payload["reaction_signature"]["schema_version"] == "1.5"
    assert payload["reaction_signature"]["topology"] == payload["reaction_topology"]


def test_same_and_different_relationship_constraints_are_enforced() -> None:
    intramolecular = featurize_reaction(INTRAMOLECULAR_CN)
    intermolecular = featurize_reaction(
        "Brc1ccccc1.CCN>>CCNc1ccccc1"
    )
    base_grammar = next(
        grammar
        for grammar in load_reaction_grammars()
        if grammar["id"] == "sp2_c_n_substitution"
    )

    same = {
        **base_grammar,
        "role_relationships": [
            {
                "roles": ["electrophile", "nucleophile"],
                "component_relation": "same",
            }
        ],
    }
    different = {
        **base_grammar,
        "role_relationships": [
            {
                "roles": ["electrophile", "nucleophile"],
                "component_relation": "different",
            }
        ],
    }

    assert len(enumerate_role_assignments(intramolecular.reactants, same)) == 1
    assert enumerate_role_assignments(intramolecular.reactants, different) == []
    assert enumerate_role_assignments(intermolecular.reactants, same) == []
    assert len(enumerate_role_assignments(intermolecular.reactants, different)) == 1
