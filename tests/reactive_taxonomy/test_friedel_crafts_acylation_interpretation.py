"""Regression coverage for the Friedel–Crafts acylation interpretation."""

from reactive_taxonomy import featurize_reaction


def test_acyl_chloride_and_benzene_reconstruct_aryl_ketone() -> None:
    result = featurize_reaction("c1ccccc1.CC(=O)Cl>>CC(=O)c1ccccc1")

    assert result.reaction_label.concise == "R–C(O)Cl + Ar–H → Ar–C(O)–R"
    assert result.reaction_label.status == "family_overlay"
    assert result.transformation_class == "aromatic_c_h_acylation"
    assert result.named_family == "friedel_crafts_acylation"
    assert result.evidence_quality == "exact_product_reconstruction"
    assert "SYMMETRY_EQUIVALENT_ASSIGNMENTS" in result.warnings


def test_two_aryl_fragments_keep_role_provenance_in_product_label() -> None:
    result = featurize_reaction("O=C(Cl)c1ccccc1.c1ccccc1>>O=C(c1ccccc1)c1ccccc1")

    assert result.reaction_label.concise == ("Ar1–C(O)Cl + Ar2–H → Ar2–C(O)–Ar1")


def test_acid_anhydride_is_an_activated_acyl_partner() -> None:
    result = featurize_reaction("CC(=O)OC(=O)C.c1ccccc1>>CC(=O)c1ccccc1")

    assert result.reaction_label.concise == ("R–C(O)OC(O)R + Ar–H → Ar–C(O)–R")
    assert result.named_family == "friedel_crafts_acylation"


def test_carboxylic_acid_is_not_promoted_to_activated_acyl_partner() -> None:
    result = featurize_reaction("CC(=O)O.c1ccccc1>>CC(=O)c1ccccc1")

    assert all(
        candidate.annotation_id != "friedel_crafts_acylation"
        for candidate in result.candidates
    )
    assert result.named_family is None


def test_non_acylated_product_does_not_verify() -> None:
    result = featurize_reaction("CC(=O)Cl.c1ccccc1>>CC(O)c1ccccc1")

    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.annotation_id == "friedel_crafts_acylation"
    ]
    assert len(candidates) == 6
    assert all(candidate.verification == "product_mismatch" for candidate in candidates)
    assert result.selected_candidate is None
    assert result.named_family is None
    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_label.status == "observed_edits"
    assert result.transformation_class == "generic_graph_transformation"


def test_product_reconstruction_resolves_aromatic_regioisomer() -> None:
    result = featurize_reaction("Cc1ccc(cc1).CC(=O)Cl>>CC(=O)c1ccc(C)cc1")

    exact = [
        candidate
        for candidate in result.candidates
        if candidate.annotation_id == "friedel_crafts_acylation"
        and candidate.verification == "exact_product_reconstruction"
    ]
    assert len(exact) == 1
    assert result.selected_candidate is not None
    assert result.named_family == "friedel_crafts_acylation"


def test_disagreeing_product_leaves_multiple_sites_unselected() -> None:
    result = featurize_reaction("Cc1ccccc1.CC(=O)Cl>>CCO")

    candidates = [
        candidate
        for candidate in result.candidates
        if candidate.annotation_id == "friedel_crafts_acylation"
    ]
    assert len(candidates) == 5
    assert result.selected_candidate is None
    assert all(candidate.verification == "product_mismatch" for candidate in candidates)


def test_mapped_acylation_agrees_with_operator() -> None:
    result = featurize_reaction(
        "[CH3:1][C:2](=[O:3])[Cl:4]."
        "[cH:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1>>"
        "[CH3:1][C:2](=[O:3])[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1"
    )

    assert result.reaction_signature is not None
    assert result.reaction_signature.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )


def test_mapping_operator_conflict_remains_visible() -> None:
    result = featurize_reaction(
        "[CH3:1][C:2](=[O:3])[Cl:4]."
        "[cH:5]1[cH:6][cH:7][c:8]([CH3:9])[cH:10][cH:11]1>>"
        "[CH3:1][C:2](=[O:3])[c:6]1[cH:5][cH:7]"
        "[c:8]([CH3:9])[cH:10][cH:11]1"
    )

    assert result.reaction_label.status == "conflicting_evidence"
    assert result.evidence_quality == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.warnings


def test_reactant_order_and_repeated_runs_are_deterministic() -> None:
    forward = featurize_reaction("c1ccccc1.CC(=O)Cl>>CC(=O)c1ccccc1")
    reversed_order = featurize_reaction("CC(=O)Cl.c1ccccc1>>CC(=O)c1ccccc1")
    repeated = featurize_reaction("c1ccccc1.CC(=O)Cl>>CC(=O)c1ccccc1")

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert repeated.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )
    assert forward.reaction_signature.signature_id == (
        repeated.reaction_signature.signature_id
    )
