"""Regression coverage for reusable high-ROI reaction grammars."""

from reactive_taxonomy import featurize_reaction


ACYL_SULFUR_COUPLING = (
    "CC(=O)S.CC(=O)SCc1cccc(I)c1.[K]>>"
    "CC(=O)SCc1cccc(SC(C)=O)c1"
)


def test_sp3_oxygen_substitution_reconstructs_product() -> None:
    result = featurize_reaction("CCCCBr.Oc1ccccc1>>CCCCOc1ccccc1")

    assert result.reaction_label.concise == "R–Br + Ar–OH → R–O–Ar"
    assert result.transformation_class == "sp3_c_o_substitution"
    assert result.named_family is None
    assert result.evidence_quality == "exact_product_reconstruction"


def test_sp3_nitrogen_and_sulfur_substitutions_reconstruct_products() -> None:
    nitrogen = featurize_reaction("CN.CCBr>>CCNC")
    sulfur = featurize_reaction("CI.SC>>CSC")

    assert nitrogen.transformation_class == "sp3_c_n_substitution"
    assert nitrogen.reaction_label.concise == "R1–Br + R2–NH2 → R1–NH–R2"
    assert sulfur.transformation_class == "sp3_c_s_substitution"
    assert sulfur.reaction_label.concise == "R1–I + R2–SH → R1–S–R2"


def test_acyl_sulfur_donor_reconstructs_sp2_c_s_coupling() -> None:
    result = featurize_reaction(ACYL_SULFUR_COUPLING)

    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.transformation_class == "sp2_c_s_substitution"
    assert result.named_family == "c_s_coupling"
    assert result.selected_candidate is not None
    assert result.selected_candidate.grammar_id == "sp2_c_s_acyl_substitution"
    assert result.reaction_label.concise == (
        "Ar–I + R–C(O)–SH → Ar–S–C(O)–R"
    )
    assert result.reaction_label is not None
    assert result.reaction_label.detailed == (
        "Ar–I + R–C(O)–SH → Ar–S–C(O)–R; edits: "
        "C–I bond cleavage; C–S bond formation; S–H loss"
    )
    sulfur_site = result.selected_candidate.role_assignments["nucleophile"]
    assert sulfur_site.availability == "deactivated"
    assert sulfur_site.canonical_signature == "XH|S|H1|C(O)R"


def test_acyl_sulfur_coupling_is_reactant_order_invariant() -> None:
    forward = featurize_reaction(ACYL_SULFUR_COUPLING)
    reversed_order = featurize_reaction(
        "CC(=O)SCc1cccc(I)c1.CC(=O)S.[K]>>"
        "CC(=O)SCc1cccc(SC(C)=O)c1"
    )

    assert forward.evidence_quality == "exact_product_reconstruction"
    assert reversed_order.evidence_quality == "exact_product_reconstruction"
    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert (
        forward.reaction_signature.signature_id
        == reversed_order.reaction_signature.signature_id
    )


def test_explicit_potassium_thioacetate_reconstructs_without_h_loss() -> None:
    result = featurize_reaction(
        "CC(=O)[S-].[K+].Ic1ccccc1>>CC(=O)Sc1ccccc1"
    )

    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.named_family == "c_s_coupling"
    assert result.selected_candidate is not None
    assert result.selected_candidate.grammar_id == "sp2_c_s_anion_substitution"
    assert result.reaction_label.concise == (
        "Ar–I + R–C(O)–S⁻ → Ar–S–C(O)–R"
    )
    assert result.reaction_label is not None
    assert "S–H loss" not in result.reaction_label.detailed
    assert [
        change.change_type
        for change in result.selected_candidate.predicted_bond_changes
    ] == ["broken", "formed"]


def test_acyl_sulfur_candidate_is_not_selected_for_wrong_product() -> None:
    result = featurize_reaction(
        "CC(=O)S.Ic1ccccc1>>c1ccccc1"
    )

    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "sp2_c_s_acyl_substitution"
    ]
    assert len(candidates) == 1
    assert candidates[0].verification == "product_mismatch"
    assert result.selected_candidate is None


def test_sp3_substitution_does_not_force_a_mechanistic_family() -> None:
    result = featurize_reaction(
        "ClC(c1ccccc1)(c1ccccc1)c1ccccc1.OCC>>CCOC(c1ccccc1)(c1ccccc1)c1ccccc1"
    )

    assert result.reaction_label.status == "exact_reconstruction"
    assert result.transformation_class == "sp3_c_o_substitution"
    assert result.named_family is None


def test_sp3_product_mismatch_is_not_selected() -> None:
    result = featurize_reaction("CCBr.CN>>CCOC")

    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "sp3_c_n_substitution"
    ]
    assert len(candidates) == 1
    assert candidates[0].verification == "product_mismatch"
    assert result.selected_candidate is None


def test_multiple_sp3_assignments_remain_unselected_when_product_disagrees() -> None:
    result = featurize_reaction("CCBr.CO.CCO>>CC")

    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "sp3_c_o_substitution"
    ]
    assert len(candidates) == 2
    assert result.selected_candidate is None
    assert all(candidate.verification == "product_mismatch" for candidate in candidates)


def test_mapped_sp3_substitution_agrees_with_operator() -> None:
    result = featurize_reaction("[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]")

    assert result.reaction_label.status == "exact_reconstruction"
    assert result.reaction_signature is not None
    assert result.reaction_signature.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )


def test_mapped_sp3_conflict_is_preserved() -> None:
    result = featurize_reaction(
        "[CH3:5][CH2:2][Br:1].[NH2:3][CH3:4]>>[CH3:2][CH2:5][NH:3][CH3:4]"
    )

    assert result.reaction_label.status == "conflicting_evidence"
    assert result.evidence_quality == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.warnings


def test_terminal_alkene_heck_reconstructs_product() -> None:
    result = featurize_reaction("Brc1ccccc1.C=C>>C=Cc1ccccc1")

    assert result.reaction_label.concise == "Ar–Br + H2C=CH2 → Ar–CH=CH2"
    assert result.transformation_class == "alkene_c_c_substitution"
    assert result.named_family == "heck"
    assert result.evidence_quality == "exact_product_reconstruction"


def test_terminal_alkene_heck_handles_unsymmetrical_partner() -> None:
    result = featurize_reaction("Brc1ccncc1.C=CC#N>>N#CC=Cc1ccncc1")

    assert result.reaction_label.concise == ("HeteroAr–Br + H2C=CHR1 → HeteroAr–CH=CHR1")
    assert result.named_family == "heck"


def test_heck_stereochemistry_is_not_invented() -> None:
    result = featurize_reaction("Ic1ccc(C=O)cc1.C=CC(=O)OC>>COC(=O)/C=C/c1ccc(C=O)cc1")

    candidate = next(
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "terminal_alkene_heck_coupling"
    )
    assert candidate.verification == "product_mismatch"
    assert result.named_family is None
    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_signature is not None
    assert any(
        change.stereo_type == "bond"
        and change.new_descriptor == "E"
        and change.change_type == "created"
        for change in result.reaction_signature.stereo_changes
    )


def test_fully_substituted_alkene_is_not_terminal_heck_partner() -> None:
    result = featurize_reaction("Brc1ccccc1.CC(C)=C(C)C>>c1ccc(C2(C)C(C)(C)C2)cc1")

    assert all(
        candidate.grammar_id != "terminal_alkene_heck_coupling"
        for candidate in result.candidates
    )


def test_mapped_heck_agrees_and_reactant_order_is_invariant() -> None:
    mapped = featurize_reaction(
        "[Br:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[CH2:8]=[CH2:9]>>"
        "[c:2]1([CH:8]=[CH2:9])[cH:3][cH:4][cH:5][cH:6][cH:7]1"
    )
    forward = featurize_reaction("Brc1ccccc1.C=C>>C=Cc1ccccc1")
    reversed_order = featurize_reaction("C=C.Brc1ccccc1>>C=Cc1ccccc1")

    assert mapped.reaction_signature is not None
    assert mapped.reaction_signature.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )
    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )


def test_activated_sp3_carbon_arylation_reconstructs_product() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.CC#N>>N#CCc1ccccc1"
    )

    assert result.transformation_class == "sp2_c_c_substitution"
    assert result.named_family is None
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.selected_candidate is not None
    assert result.selected_candidate.grammar_id == (
        "sp2_c_activated_c_substitution"
    )
    assert result.product_connection is not None
    assert result.product_connection.connection_type == "C_C"


def test_activated_carbon_grammar_rejects_unactivated_alkane() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.CC>>CCc1ccccc1"
    )

    assert all(
        candidate.grammar_id != "sp2_c_activated_c_substitution"
        for candidate in result.candidates
    )


def test_aromatic_ch_arylation_reconstructs_biaryl_product() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert result.transformation_class == "sp2_c_c_substitution"
    assert result.named_family is None
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.selected_candidate is not None
    assert result.selected_candidate.grammar_id == (
        "sp2_c_aromatic_ch_substitution"
    )
    assert result.reaction_label.concise == "Ar1–Br + Ar2–H → Ar1–Ar2"
    assert result.product_connection is not None
    assert result.product_connection.connection_type == "C_C"


def test_aromatic_ch_arylation_is_reactant_order_invariant() -> None:
    forward = featurize_reaction(
        "Brc1ccccc1.c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    reversed_order = featurize_reaction(
        "c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )


def test_aromatic_ch_arylation_does_not_select_product_mismatch() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.c1ccccc1>>c1ccccc1"
    )

    direct_arylation = [
        candidate
        for candidate in result.candidates
        if candidate.grammar_id == "sp2_c_aromatic_ch_substitution"
    ]
    assert direct_arylation
    assert all(
        candidate.verification != "exact_product_reconstruction"
        for candidate in direct_arylation
    )
    assert result.selected_candidate is None
