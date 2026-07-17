"""Regression tests for conservative unmapped scaffold correspondence."""

from reactive_taxonomy import featurize_reaction


def test_unmapped_alkene_hydrogenation_receives_inferred_label() -> None:
    result = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )

    assert result.valid
    assert result.evidence_quality == "unique_scaffold_correspondence"
    assert result.reaction_label == "Ar–CH=CH2 → Ar–CH2–CH3"
    assert result.reaction_label_status == "generic_pattern"
    assert result.display_label is not None
    assert result.display_label.pattern_id == "hydrogenation"
    assert result.display_label.transformation_label == "C=C hydrogenation"
    assert result.display_label.reactant_context_label == "Ar–CH=CH2"
    assert result.display_label.product_context_label == "Ar–CH2–CH3"
    assert result.display_label.confidence == 0.85
    assert result.reaction_signature is not None
    assert result.reaction_signature.product_transformation is not None
    assert (
        result.reaction_signature.product_transformation.concise_label
        == "Ar–CH2–CH3"
    )
    assert "INFERRED_ATOM_CORRESPONDENCE" in result.warnings


def test_unmapped_carbonyl_redox_receives_generic_labels() -> None:
    reduction = featurize_reaction(
        "c1ccc(C=O)cc1.[H][H]>>c1ccc(CO)cc1"
    )
    oxidation = featurize_reaction(
        "CC(O)c1ccccc1.[O]>>CC(=O)c1ccccc1"
    )

    assert reduction.reaction_label == "Ar–CH=O → Ar–CH2–OH"
    assert reduction.display_label is not None
    assert reduction.display_label.pattern_id == "heteroatom_bond_reduction"
    assert reduction.display_label.transformation_label == "C=O reduction"
    assert oxidation.reaction_label == "Ar–CH(R)–OH → Ar–C(R)=O"
    assert oxidation.display_label is not None
    assert oxidation.display_label.pattern_id == "heteroatom_bond_oxidation"


def test_unmapped_complete_alkyne_hydrogenation_receives_generic_label() -> None:
    result = featurize_reaction(
        "C#Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )

    assert result.reaction_label == "Ar–C≡CH → Ar–CH2–CH3"
    assert result.display_label is not None
    assert result.display_label.pattern_id == "complete_alkyne_hydrogenation"
    assert result.display_label.structural_label == (
        "C≡C → C–C; 4 × H gain at C"
    )


def test_symmetry_equivalent_correspondence_is_order_invariant() -> None:
    forward = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )
    reversed_order = featurize_reaction(
        "[H][H].C=Cc1ccccc1>>CCc1ccccc1"
    )

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )
    assert forward.reaction_label == reversed_order.reaction_label


def test_contextual_label_style_changes_display_only() -> None:
    unicode_result = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )
    ascii_result = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1", label_style="ascii"
    )

    assert unicode_result.reaction_label == "Ar–CH=CH2 → Ar–CH2–CH3"
    assert ascii_result.reaction_label == "Ar-CH=CH2 -> Ar-CH2-CH3"
    assert unicode_result.reaction_signature is not None
    assert ascii_result.reaction_signature is not None
    assert unicode_result.reaction_signature.signature_id == (
        ascii_result.reaction_signature.signature_id
    )


def test_multiple_order_changes_keep_literal_edit_fallback() -> None:
    result = featurize_reaction(
        "[CH2:1]=[CH:2][CH:3]=[CH2:4]>>"
        "[CH3:1][CH2:2][CH2:3][CH3:4]"
    )

    assert result.display_label is not None
    assert result.display_label.contextual_label is None
    assert result.display_label.pattern_id is None
    assert result.reaction_label == "2 × C=C → C–C; 4 × H gain at C"


def test_chemically_distinct_correspondences_remain_ambiguous() -> None:
    result = featurize_reaction("CCCBr>>CC=C")

    assert result.valid
    assert result.evidence_quality == "ambiguous_atom_correspondence"
    assert result.reaction_label is None
    assert result.reaction_signature is None
    assert "AMBIGUOUS_SCAFFOLD_CORRESPONDENCE:2" in result.warnings


def test_multi_substrate_assembly_is_outside_scaffold_fallback() -> None:
    result = featurize_reaction("CC.CN>>CCN")

    assert result.valid
    assert result.evidence_quality == "unresolved"
    assert result.reaction_label is None
    assert result.reaction_signature is None
    assert "SCAFFOLD_CORRESPONDENCE_REQUIRES_ONE_SUBSTRATE" in result.warnings
