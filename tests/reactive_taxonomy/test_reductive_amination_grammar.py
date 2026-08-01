"""Regression coverage for carbonyl–amine reductive coupling grammar."""

from reactive_taxonomy import featurize_reaction


def test_aldehyde_primary_amine_exact_reconstruction() -> None:
    result = featurize_reaction("O=Cc1ccncc1.Nc1ccccc1>>c1ccc(NCc2ccncc2)cc1")

    assert result.valid
    assert result.reaction_label.concise == ("HeteroAr–CH=O + Ar–NH2 → HeteroAr–CH2–NH–Ar")
    assert result.reaction_label.status == "family_overlay"
    assert result.transformation_class == "carbonyl_c_n_reductive_coupling"
    assert result.named_family == "reductive_amination"
    assert result.selected_candidate is not None
    assert result.selected_candidate.grammar_id == ("carbonyl_amine_reductive_coupling")
    assert result.reaction_signature is not None
    assert result.reaction_signature.product_transformation is not None
    assert result.reaction_signature.product_transformation.concise_label is None


def test_ketone_primary_amine_exact_reconstruction() -> None:
    result = featurize_reaction("CC(=O)c1ccccc1.Nc1ccccc1>>CC(Nc1ccccc1)c1ccccc1")

    assert result.reaction_label.concise == ("Ar1–C(R)=O + Ar2–NH2 → Ar1–CH(R)–NH–Ar2")
    assert result.evidence_quality == "exact_product_reconstruction"


def test_formaldehyde_primary_amine_exact_reconstruction() -> None:
    result = featurize_reaction("C=O.Nc1ccccc1>>CNc1ccccc1")

    assert result.reaction_label.concise == "H2C=O + Ar–NH2 → H3C–NH–Ar"
    assert result.named_family == "reductive_amination"


def test_aldehyde_secondary_amine_exact_reconstruction() -> None:
    result = featurize_reaction("O=Cc1ccccc1.CNC>>CN(C)Cc1ccccc1")

    assert result.reaction_label.concise == "Ar–CH=O + R1R2–NH → Ar–CH2–NR1R2"
    assert result.named_family == "reductive_amination"


def test_tertiary_amine_is_not_a_reductive_amination_partner() -> None:
    result = featurize_reaction("O=Cc1ccccc1.CN(C)C>>C[N+](C)(C)Cc1ccccc1")

    assert all(
        candidate.grammar_id != "carbonyl_amine_reductive_coupling"
        for candidate in result.candidates
    )
    assert result.named_family is None


def test_amide_nitrogen_is_not_a_reductive_amination_partner() -> None:
    result = featurize_reaction("O=Cc1ccccc1.NC(=O)c1ccccc1>>O=C(NCc1ccccc1)c1ccccc1")

    assert all(
        candidate.grammar_id != "carbonyl_amine_reductive_coupling"
        for candidate in result.candidates
    )


def test_imine_product_does_not_verify_as_reductive_amination() -> None:
    result = featurize_reaction("O=Cc1ccncc1.Nc1ccccc1>>c1ccc(N=Cc2ccncc2)cc1")

    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "carbonyl_amine_reductive_coupling"
    ]
    assert len(candidates) == 1
    assert candidates[0].verification == "product_mismatch"
    assert result.named_family is None
    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_label.status == "generic_pattern"
    assert result.transformation_class == "substitution"


def test_multiple_carbonyl_partners_remain_unselected_without_verification() -> None:
    result = featurize_reaction("O=Cc1ccccc1.CC=O.Nc1ccccc1>>CC")

    reductive_candidates = [
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "carbonyl_amine_reductive_coupling"
    ]
    assert len(reductive_candidates) == 2
    assert result.selected_candidate is None
    assert result.reaction_label.status == "product_contradicted_reactants"


def test_mapped_reductive_amination_agrees_with_operator() -> None:
    result = featurize_reaction(
        "[O:1]=[CH:2]c1ccncc1.[NH2:3]c1ccccc1>>[CH2:2]([NH:3]c1ccccc1)c1ccncc1"
    )

    assert result.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )
    assert result.reaction_signature is not None
    assert result.reaction_signature.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )


def test_mapping_operator_conflict_is_visible() -> None:
    result = featurize_reaction(
        "[O:1]=[CH:2][c:4]1ccccc1.[NH2:3]c1ccccc1>>[CH2:4]([NH:3]c1ccccc1)[c:2]1ccccc1"
    )

    assert result.reaction_label.status == "conflicting_evidence"
    assert result.evidence_quality == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.warnings


def test_reactant_order_and_repeated_runs_are_deterministic() -> None:
    forward = featurize_reaction("O=Cc1ccncc1.Nc1ccccc1>>c1ccc(NCc2ccncc2)cc1")
    reversed_order = featurize_reaction("Nc1ccccc1.O=Cc1ccncc1>>c1ccc(NCc2ccncc2)cc1")
    repeated = featurize_reaction("O=Cc1ccncc1.Nc1ccccc1>>c1ccc(NCc2ccncc2)cc1")

    assert forward.reaction_label == reversed_order.reaction_label
    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert repeated.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )
    assert forward.reaction_signature.signature_id == (
        repeated.reaction_signature.signature_id
    )
