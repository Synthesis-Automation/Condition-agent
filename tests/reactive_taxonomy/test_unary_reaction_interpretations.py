"""Regression coverage for high-ROI unary transformation interpretations."""

import pytest

from reactive_taxonomy import featurize_reaction


@pytest.mark.parametrize(
    (
        "reaction",
        "annotation_id",
        "transformation_class",
        "order_change",
        "hydrogen_change_count",
    ),
    (
        (
            "CC(=O)C>>CC(O)C",
            "carbonyl_reduction",
            "carbonyl_reduction",
            "C-O:DOUBLE>SINGLE",
            2,
        ),
        (
            "CC(O)C>>CC(=O)C",
            "alcohol_oxidation",
            "alcohol_oxidation",
            "C-O:SINGLE>DOUBLE",
            2,
        ),
        (
            "CC=CC>>CCCC",
            "alkene_hydrogenation",
            "alkene_reduction",
            "C-C:DOUBLE>SINGLE",
            2,
        ),
        (
            "CC#CC>>CC=CC",
            "alkyne_partial_hydrogenation",
            "alkyne_partial_reduction",
            "C-C:TRIPLE>DOUBLE",
            2,
        ),
        (
            "CC#CC>>CCCC",
            "alkyne_hydrogenation",
            "alkyne_reduction",
            "C-C:TRIPLE>SINGLE",
            4,
        ),
    ),
)
def test_unary_bond_order_interpretation_reconstructs_exact_product(
    reaction: str,
    annotation_id: str,
    transformation_class: str,
    order_change: str,
    hydrogen_change_count: int,
) -> None:
    result = featurize_reaction(reaction)

    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.selected_candidate is not None
    assert result.selected_candidate.annotation_id == annotation_id
    assert result.transformation_class == transformation_class
    assert result.named_family is None
    assert result.reaction_signature is not None
    assert result.reaction_signature.order_changes == (order_change,)
    assert len(result.reaction_signature.hydrogen_changes) == (
        hydrogen_change_count
    )
    assert result.reaction_signature.product_transformation is not None
    assert result.reaction_signature.product_transformation.exact_product_verified


def test_sufex_sulfonate_reconstructs_without_source_label_routing() -> None:
    result = featurize_reaction(
        "Oc1ccccc1.O=S(=O)(F)CC>>O=S(=O)(CC)Oc1ccccc1"
    )

    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.selected_candidate is not None
    assert result.selected_candidate.annotation_id == "sulfonate_formation"
    assert result.transformation_class == "sulfonate_formation"
    assert result.named_family == "sulfonylation"
    assert result.reaction_label.concise == (
        "R–S(O)2F + Ar–OH → R–S(O)2–O–Ar"
    )
    assert result.reaction_signature is not None
    assert result.reaction_signature.formed_bond_types == ("O-S:SINGLE",)
    assert result.reaction_signature.broken_bond_types == ("F-S:SINGLE",)


def test_unary_operator_does_not_select_a_mismatched_product() -> None:
    result = featurize_reaction("CC(=O)C>>CC")
    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.annotation_id == "carbonyl_reduction"
    ]

    assert len(candidates) == 1
    assert candidates[0].verification == "product_mismatch"
    assert result.selected_candidate is None


def test_tertiary_alcohol_is_not_invented_as_a_carbonyl() -> None:
    result = featurize_reaction("CC(C)(C)O>>CC(C)=O")
    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.annotation_id == "alcohol_oxidation"
    ]

    assert len(candidates) == 1
    assert candidates[0].verification == "construction_failed"
    assert result.selected_candidate is None


def test_equivalent_carbonyl_sites_are_deterministic() -> None:
    result = featurize_reaction(
        "CC(=O)CC(C)=O>>CC(O)CC(C)=O"
    )

    assert result.selected_candidate is not None
    assert result.selected_candidate.annotation_id == "carbonyl_reduction"
    assert result.evidence_quality == "exact_product_reconstruction"
    assert "SYMMETRY_EQUIVALENT_ASSIGNMENTS" in result.warnings


def test_mapped_carbonyl_reduction_agrees_with_operator() -> None:
    result = featurize_reaction(
        "[CH3:1][C:2](=[O:3])[CH3:4]>>"
        "[CH3:1][CH:2]([OH:3])[CH3:4]"
    )

    assert result.reaction_signature is not None
    assert result.reaction_signature.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )


def test_mapping_operator_conflict_remains_review_evidence() -> None:
    result = featurize_reaction(
        "[CH3:1][C:2](=[O:3])[CH3:4]>>"
        "[CH3:2][CH:1]([OH:3])[CH3:4]"
    )

    assert result.reaction_signature is not None
    assert result.evidence_quality == "conflicting_edit_evidence"
    assert result.reaction_label.status == "conflicting_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.warnings
