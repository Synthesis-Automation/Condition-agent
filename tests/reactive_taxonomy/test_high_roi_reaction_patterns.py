"""High-ROI pharmaceutical reaction-pattern regressions."""

from __future__ import annotations

import json
from pathlib import Path

from reactive_taxonomy import featurize_reaction, identify_reaction_patterns


def _patterns(reaction: str) -> dict[str, object]:
    result = featurize_reaction(reaction)
    assert result.interpretation is not None
    return {
        match.pattern_id: match for match in result.interpretation.pattern_matches
    }


def test_sp2_c_n_pattern_retains_named_family_ambiguity() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.CN>>CNc1ccccc1"
    )

    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == "sp2_c_n_substitution_like"
    assert result.named_family is None
    assert result.compatible_named_families == (
        "buchwald_hartwig_c_n",
        "snar_amination",
        "ullmann_c_n",
    )
    match = result.interpretation.pattern_matches[0]
    assert match.requires_condition_evidence
    assert result.interpretation.co_occurring_pattern_ids == ()
    assert result.reaction_label is not None
    assert result.reaction_label.concise == "C(sp²)–N bond formation"
    assert "named family requires condition evidence" in result.reaction_label.detailed


def test_multi_site_sp2_c_n_pattern_preserves_both_events() -> None:
    result = featurize_reaction(
        "Brc1ccc(Br)cc1.CC(N)=O"
        ">>CC(=O)Nc1ccc(NC(C)=O)cc1"
    )

    assert result.reaction_completeness is not None
    assert result.reaction_completeness.status == "verified"
    assert (
        result.reaction_completeness.suspected_insufficient_reactant_multiplicity
    )
    assert result.reactants[2].inferred_copy_of_component_index == 1
    assert result.reaction_signature is not None
    assert result.reaction_signature.event_count == 2
    assert result.reaction_signature.formed_bond_types == (
        "C-N:SINGLE",
        "C-N:SINGLE",
    )
    assert result.reaction_signature.broken_bond_types == (
        "Br-C:SINGLE",
        "Br-C:SINGLE",
    )
    assert result.interpretation is not None
    pattern = next(
        match
        for match in result.interpretation.pattern_matches
        if match.pattern_id == "sp2_c_n_substitution_like"
    )
    assert pattern.occurrence_count == 2
    assert pattern.matched_edit_indices == (0, 1, 2, 3, 4, 5)
    assert result.reaction_core is not None
    assert result.reaction_core.generic_label == (
        "2 × Ar–Br + 2 × H–N → 2 × Ar–N"
    )
    assert result.reaction_label is not None
    assert result.reaction_label.concise == (
        "2 × (C–Br + N–H → C–N)"
    )


def test_explicit_multi_site_reagent_copies_avoid_inference_warning() -> None:
    result = featurize_reaction(
        "Brc1ccc(Br)cc1.CC(N)=O.CC(N)=O"
        ">>CC(=O)Nc1ccc(NC(C)=O)cc1"
    )

    assert result.reaction_signature is not None
    assert result.reaction_signature.event_count == 2
    assert "INFERRED_REACTANT_MULTIPLICITY" not in result.warnings
    assert all(
        component.inferred_copy_of_component_index is None
        for component in result.reactants
    )

    shorthand = featurize_reaction(
        "Brc1ccc(Br)cc1.CC(N)=O"
        ">>CC(=O)Nc1ccc(NC(C)=O)cc1"
    )
    assert shorthand.reaction_signature is not None
    assert (
        shorthand.reaction_signature.signature_id
        == result.reaction_signature.signature_id
    )


def test_tandem_c_n_coupling_and_deprotection_are_both_retained() -> None:
    result = featurize_reaction(
        "CC(C)(C)OC(=O)NCCN."
        "O=c1oc2cc(Br)ccc2cc1-c1ccccc1"
        ">>NCCNc1ccc2cc(-c3ccccc3)c(=O)oc2c1"
    )

    assert result.reaction_signature is not None
    assert "INFERRED_SINGLE_CUT_FRAGMENT_CORRESPONDENCE" in result.warnings
    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == (
        "sp2_c_n_substitution_like"
    )
    assert result.interpretation.co_occurring_pattern_ids == (
        "sp2_c_n_substitution_like",
        "amine_deprotection_like",
    )
    assert result.reaction_label is not None
    assert result.reaction_label.concise == (
        "C(sp²)–N bond formation + amine deprotection"
    )
    assert "temporal order is not inferred" in result.reaction_label.detailed


def test_public_pattern_identification_accepts_reaction_smiles() -> None:
    interpretation = identify_reaction_patterns(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert interpretation is not None
    assert interpretation.primary_pattern_id == "organoboron_c_c_coupling_like"
    assert interpretation.named_family == "suzuki_miyaura"


def test_sp2_c_n_pattern_excludes_acylation_and_reductive_amination() -> None:
    amide = _patterns("CC(=O)Cl.CN>>CC(=O)NC")
    reductive = _patterns("CC=O.CN>>CCNC")

    assert "amide_formation_like" in amide
    assert "sp2_c_n_substitution_like" not in amide
    assert "reductive_amination_like" in reductive
    assert "sp2_c_n_substitution_like" not in reductive


def test_sp3_substitution_distinguishes_n_and_o_installation() -> None:
    nitrogen = _patterns("CCBr.CN>>CCNC")
    oxygen = _patterns("CCBr.CO>>CCOC")

    assert "sp3_c_n_substitution_like" in nitrogen
    assert "sp3_c_o_substitution_like" not in nitrogen
    assert "sp3_c_o_substitution_like" in oxygen
    assert "sp3_c_n_substitution_like" not in oxygen


def test_chan_lam_like_connectivity_is_distinct_from_suzuki_and_sp2_substitution() -> None:
    chan_lam = featurize_reaction("OB(O)c1ccccc1.CN>>CNc1ccccc1")
    suzuki = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert chan_lam.interpretation is not None
    assert chan_lam.interpretation.primary_pattern_id == (
        "organoboron_c_n_coupling_like"
    )
    assert chan_lam.named_family is None
    assert chan_lam.compatible_named_families == ("chan_lam_c_n",)
    assert "sp2_c_n_substitution_like" not in {
        match.pattern_id for match in chan_lam.interpretation.pattern_matches
    }
    assert chan_lam.reaction_label is not None
    assert chan_lam.reaction_label.concise == "organoboron C–N coupling"

    assert suzuki.interpretation is not None
    assert suzuki.interpretation.primary_pattern_id == (
        "organoboron_c_c_coupling_like"
    )
    assert suzuki.named_family == "suzuki_miyaura"


def test_carbonyl_alpha_c_h_coupling_is_not_mislabeled_as_heck() -> None:
    alpha_c_h = featurize_reaction(
        "CC(=O)c1ccc(Br)cc1.O=C1CCCC1"
        ">>CC(=O)c1ccc(C2CCCC2=O)cc1"
    )
    heck = featurize_reaction("Brc1ccccc1.C=C>>C=Cc1ccccc1")

    assert alpha_c_h.interpretation is not None
    assert alpha_c_h.interpretation.primary_pattern_id == (
        "carbonyl_alpha_c_h_sp2_c_c_coupling_like"
    )
    assert alpha_c_h.named_family == "carbonyl_alpha_c_h_coupling"
    assert "heck_coupling_like" not in {
        match.pattern_id for match in alpha_c_h.interpretation.pattern_matches
    }
    assert alpha_c_h.reaction_label is not None
    assert alpha_c_h.reaction_label.concise == (
        "C(sp²)–C coupling at carbonyl α-C–H"
    )

    assert heck.interpretation is not None
    assert heck.interpretation.primary_pattern_id == "heck_coupling_like"
    assert heck.named_family == "heck"
    assert "carbonyl_alpha_c_h_sp2_c_c_coupling_like" not in {
        match.pattern_id for match in heck.interpretation.pattern_matches
    }


def test_acyl_and_sulfonyl_products_receive_specific_patterns() -> None:
    examples = {
        "CC(=O)Cl.CN>>CC(=O)NC": "amide_formation_like",
        "CC(=O)Cl.CO>>CC(=O)OC": "ester_formation_like",
        "COC(=O)Cl.CN>>COC(=O)NC": "carbamate_formation_like",
        "COC(=O)Cl.CO>>COC(=O)OC": "carbonate_formation_like",
        "CN=C=O.CN>>CNC(=O)NC": "urea_formation_like",
        "CS(=O)(=O)Cl.CN>>CS(=O)(=O)NC": "sulfonamide_formation_like",
        "CS(=O)(=O)Cl.CO>>CS(=O)(=O)OC": "o_sulfonylation_like",
    }

    for reaction, expected in examples.items():
        assert expected in _patterns(reaction)


def test_common_redox_patterns_are_element_and_direction_specific() -> None:
    assert "alkene_hydrogenation_like" in _patterns("C=C>>CC")
    assert "carbonyl_reduction_like" in _patterns("CC(=O)C>>CC(O)C")
    assert "alcohol_oxidation_like" in _patterns("CCO>>CC=O")
    assert "nitro_reduction_like" in _patterns(
        "O=[N+]([O-])c1ccccc1>>Nc1ccccc1"
    )
    assert "carbonyl_reduction_like" not in _patterns("C=C>>CC")
    assert "alcohol_oxidation_like" not in _patterns("CC(=O)C>>CC(O)C")


def test_deprotection_hydrolysis_and_fgi_patterns_use_observed_edits() -> None:
    examples = {
        (
            "[CH3:1][NH:2][C:3](=[O:4])[O:5][C:6]([CH3:7])"
            "([CH3:8])[CH3:9]>>[CH3:1][NH2:2]"
        ): "amine_deprotection_like",
        (
            "[CH3:1][O:2][Si:3]([CH3:4])([CH3:5])[CH3:6]"
            ">>[CH3:1][OH:2]"
        ): "alcohol_deprotection_like",
        (
            "[CH3:1][C:2](=[O:3])[O:4][CH3:5]"
            ">>[CH3:1][C:2](=[O:3])[OH:4]"
        ): "ester_hydrolysis_like",
        (
            "[CH3:1][C:2](=[O:3])[OH:4].[Cl:5]"
            ">>[CH3:1][C:2](=[O:3])[Cl:5]"
        ): "acid_chloride_formation_like",
        (
            "[CH3:1][CH2:2][OH:3].[Br:4]"
            ">>[CH3:1][CH2:2][Br:4]"
        ): "alcohol_to_halide_like",
        (
            "[c:1]1([Br:2])[cH:3][cH:4][cH:5][cH:6][cH:7]1."
            "[C-:8]#[N:9]>>[c:1]1([C:8]#[N:9])[cH:3][cH:4]"
            "[cH:5][cH:6][cH:7]1"
        ): "cyanation_like",
        (
            "[CH3:1][S:2][CH3:3].[O:4]"
            ">>[CH3:1][S:2](=[O:4])[CH3:3]"
        ): "sulfide_oxidation_like",
    }

    for reaction, expected in examples.items():
        assert expected in _patterns(reaction)


def test_pharma_catalog_is_complete_and_evaluation_only() -> None:
    path = (
        Path(__file__).parents[2]
        / "condition_recommender"
        / "definitions"
        / "pharma_pattern_catalog.v1.json"
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    classes = payload["classes"]

    assert [item["rank"] for item in classes] == list(range(1, 45))
    assert len({item["id"] for item in classes}) == 44
    assert "never an interpretation input" in payload["purpose"]
    assert next(item for item in classes if item["rank"] == 31)[
        "implementation"
    ] == "excluded_operation"
