"""Regression coverage for general alcohol C–O displacement rewrites."""

import pytest

from reactive_taxonomy import featurize_reaction


MITSUNOBU_LIKE_ETHERIFICATION = (
    "CNCC[C@H](O)c1ccccc1.Cc1ccccc1O"
    ">>CNCC[C@@H](Oc1ccccc1C)c1ccccc1"
)


def test_inverting_alcohol_etherification_is_exactly_reconstructed() -> None:
    analysis = featurize_reaction(MITSUNOBU_LIKE_ETHERIFICATION)

    assert analysis.evidence_quality == "exact_product_reconstruction"
    assert analysis.selected_candidate is not None
    assert analysis.selected_candidate.grammar_id == "alcohol_c_o_displacement"
    assert analysis.selected_candidate.rewrite_outcome_id == "inversion_if_defined"
    assert analysis.transformation_class == "sp3_c_o_substitution"
    assert analysis.named_family is None
    assert analysis.reaction_completeness is not None
    assert analysis.reaction_completeness.status == "verified"
    assert analysis.reaction_signature is not None
    assert {
        edit.edit_type for edit in analysis.reaction_signature.edits
    } == {"broken", "formed", "hydrogen_change"}
    assert [
        (
            change.old_descriptor,
            change.new_descriptor,
            change.change_type,
        )
        for change in analysis.reaction_signature.stereo_changes
    ] == [("S", "R", "inverted")]
    assert (
        analysis.reaction_signature.product_transformation is not None
    )
    assert (
        analysis.reaction_signature.product_transformation.stereo_changes
        == analysis.reaction_signature.stereo_changes
    )


def test_observed_retention_is_not_forced_into_an_inversion_outcome() -> None:
    reaction = (
        "CNCC[C@H](O)c1ccccc1.Cc1ccccc1O"
        ">>CNCC[C@H](Oc1ccccc1C)c1ccccc1"
    )

    analysis = featurize_reaction(reaction)

    assert analysis.evidence_quality == "exact_product_reconstruction"
    assert analysis.selected_candidate is not None
    assert analysis.selected_candidate.rewrite_outcome_id == "retention_if_defined"
    assert [
        change.change_type
        for change in analysis.reaction_signature.stereo_changes
    ] == ["retained"]


@pytest.mark.parametrize(
    ("reaction", "grammar_id", "connection_type"),
    [
        ("CC(O)C.CN>>CNC(C)C", "alcohol_c_n_displacement", "C_N"),
        ("CC(O)C.CO>>COC(C)C", "alcohol_c_o_displacement", "C_O"),
        ("CC(O)C.CS>>CSC(C)C", "alcohol_c_s_displacement", "C_S"),
    ],
)
def test_alcohol_displacement_is_shared_across_xh_partner_elements(
    reaction: str,
    grammar_id: str,
    connection_type: str,
) -> None:
    analysis = featurize_reaction(reaction)

    assert analysis.evidence_quality == "exact_product_reconstruction"
    assert analysis.selected_candidate is not None
    assert analysis.selected_candidate.grammar_id == grammar_id
    assert analysis.product_connection is not None
    assert analysis.product_connection.connection_type == connection_type


def test_achiral_tertiary_alcohol_displacement_remains_structurally_generic() -> None:
    analysis = featurize_reaction("CC(C)(C)O.CO>>COC(C)(C)C")

    assert analysis.evidence_quality == "exact_product_reconstruction"
    assert analysis.selected_candidate is not None
    assert analysis.selected_candidate.grammar_id == "alcohol_c_o_displacement"
    assert analysis.reaction_signature is not None
    assert analysis.reaction_signature.stereo_changes == ()


def test_reactant_order_does_not_change_the_displacement_signature() -> None:
    reversed_reactants = (
        "Cc1ccccc1O.CNCC[C@H](O)c1ccccc1"
        ">>CNCC[C@@H](Oc1ccccc1C)c1ccccc1"
    )

    forward = featurize_reaction(MITSUNOBU_LIKE_ETHERIFICATION)
    reversed_analysis = featurize_reaction(reversed_reactants)

    assert forward.reaction_signature is not None
    assert reversed_analysis.reaction_signature is not None
    assert (
        forward.reaction_signature.signature_id
        == reversed_analysis.reaction_signature.signature_id
    )


def test_supplied_mapping_can_contradict_the_predicted_oxygen_origin() -> None:
    reaction = (
        "[CH3:1][CH:2]([OH:3])[CH3:4].[OH:5][CH3:6]"
        ">>[CH3:1][CH:2]([O:3][CH3:6])[CH3:4]"
    )

    analysis = featurize_reaction(reaction)

    assert analysis.selected_candidate is not None
    assert analysis.evidence_quality == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in analysis.warnings
    assert analysis.reaction_signature is not None
    assert {
        (
            edit.edit_type,
            edit.atom_1.atom_map_number,
            edit.atom_2.atom_map_number if edit.atom_2 else None,
        )
        for edit in analysis.reaction_signature.edits
    } >= {
        ("formed", 3, 6),
        ("broken", 5, 6),
    }
