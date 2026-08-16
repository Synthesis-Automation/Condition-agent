"""Forward library, product prediction, and competition regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import build_generic_library
from forward_synthesis import (
    ForwardOperatorLibrary,
    assess_proposed_step,
    build_forward_library,
    load_forward_library,
    predict_products,
    save_forward_library,
)


S_ALKYLATION = "Cc1[nH]cnc1CCl.NCCS>>NCCSCc1c(C)[nH]cn1"
N_ALKYLATION = "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1"


def _library(*reactions: str) -> ForwardOperatorLibrary:
    generic = build_generic_library(
        (
            {
                "reaction_id": f"reaction-{index}",
                "reference_id": f"reference-{index}",
                "reaction_smiles": reaction,
            }
            for index, reaction in enumerate(reactions)
        ),
        levels=("L1", "L2"),
        admission_mode="data_driven",
    )
    return build_forward_library(generic)


def test_forward_library_requires_independent_source_round_trip() -> None:
    library = _library("CCBr.N>>CCN")

    assert library.source_template_count == 2
    assert library.admitted_operator_count == 2
    assert library.rejection_counts == {}
    assert all(operator.forward_smarts for operator in library.operators)


def test_forward_library_round_trips_through_json_gzip(tmp_path) -> None:
    library = _library("CCBr.N>>CCN")
    path = tmp_path / "forward-library.json.gz"

    save_forward_library(library, path)
    restored = load_forward_library(path)

    assert restored == library


def test_blind_prediction_recovers_product_without_target_input() -> None:
    library = _library("CCBr.N>>CCN")

    result = predict_products("N.CCBr", library)

    assert result.valid
    assert result.status == "predicted"
    assert tuple(candidate.product_smiles for candidate in result.candidates) == (
        "CCN",
    )
    assert result.candidates[0].operator_edit_agreement
    assert result.candidates[0].reverse_round_trip
    assert len(result.candidates[0].alternative_pathway_ids) == 1
    assert result.diagnostics.valid_pathway_count == 2
    assert "CONDITIONS_NOT_SUPPLIED_PRODUCTS_ARE_POSSIBILITIES" in result.warnings


def test_blind_prediction_reports_a_truncated_operator_search() -> None:
    library = _library(S_ALKYLATION, N_ALKYLATION)

    result = predict_products(
        "Cc1[nH]cnc1CCl.NCCS",
        library,
        max_operators_to_apply=1,
    )

    assert result.diagnostics.indexed_operator_count > 1
    assert "FORWARD_OPERATOR_SEARCH_TRUNCATED" in result.warnings


def test_competing_endpoint_products_are_retained_and_route_is_competitive() -> None:
    library = _library(S_ALKYLATION, N_ALKYLATION)
    starting_materials = "Cc1[nH]cnc1CCl.NCCS"

    prediction = predict_products(starting_materials, library)
    assessment = assess_proposed_step(
        starting_materials,
        "NCCSCc1c(C)[nH]cn1",
        library,
    )

    assert prediction.diagnostics.unique_product_count == 2
    assert len(prediction.candidates) == 2
    assert {group.competition_level for group in prediction.competition_groups} == {
        "operator",
        "site",
        "product",
    }
    assert assessment.intended_product_rank == 2
    assert assessment.disposition == "competitive"
    assert assessment.validity == "structurally_supported_with_competition"
    assert {
        check.check_id: check.status for check in assessment.checks
    }["competing_product_risk"] == "warning"
    assert assessment.best_competitor_product is not None


def test_missing_operator_hint_is_out_of_scope_not_a_chemical_contradiction() -> None:
    library = _library("CCBr.N>>CCN")

    assessment = assess_proposed_step(
        "CCBr.N",
        "CCN",
        library,
        operator_hint="OP1:not-in-forward-library",
    )

    assert assessment.targeted_replay_status == "operator_hint_not_found"
    assert assessment.disposition == "out_of_scope"
    assert assessment.validity == "out_of_scope"
    assert "RETROSYNTHESIS_OPERATOR_NOT_FORWARD_ADMITTED" in assessment.warnings
    checks = {check.check_id: check.status for check in assessment.checks}
    assert checks["targeted_operator_replay"] == "not_evaluated"
    assert checks["blind_target_recovery"] == "pass"


def test_targeted_replay_does_not_rescue_blind_operator_mismatch() -> None:
    library = _library("CCBr.N>>CCN")
    corrupted = replace(
        library.operators[0],
        edit_tokens=("formed:C-O:NONE>SINGLE",),
    )
    bad_library = replace(
        library,
        operators=(corrupted,),
        admitted_operator_count=1,
        precursor_index=replace(
            library.precursor_index,
            component_count_to_operator_ids={2: (corrupted.forward_operator_id,)},
            operator_required_atomic_numbers={
                corrupted.forward_operator_id: library.precursor_index.operator_required_atomic_numbers[
                    library.operators[0].forward_operator_id
                ]
            },
        ),
    )

    result = predict_products("CCBr.N", bad_library)

    assert not result.candidates
    assert result.diagnostics.operator_edit_mismatch_count == 1


def test_invalid_query_returns_typed_failure() -> None:
    result = predict_products("not-smiles", _library("CCBr.N>>CCN"))

    assert not result.valid
    assert result.error == "INVALID_STARTING_MATERIALS"


def test_recipe_compatibility_is_audited_not_presented_as_yield() -> None:
    result = predict_products(
        "CCBr.N",
        _library("CCBr.N>>CCN"),
        recipe={"temperature_c": 25.0, "atmosphere": "nitrogen"},
    )

    assert result.conditions_supplied
    assert result.candidates[0].recipe_evidence.evaluated
    assert result.candidates[0].recipe_evidence.compatible
    assert result.candidates[0].recipe_evidence.score is not None


def test_bromoaniline_self_coupling_is_an_explicit_competing_pathway() -> None:
    library = _library(
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    )

    included = predict_products("Brc1ccc(N)cc1", library)
    excluded = predict_products(
        "Brc1ccc(N)cc1",
        library,
        include_self_reactions=False,
    )

    assert tuple(item.product_smiles for item in included.candidates) == (
        "Nc1ccc(Nc2ccc(Br)cc2)cc1",
    )
    candidate = included.candidates[0]
    assert candidate.uses_virtual_copies is True
    assert candidate.reactant_stoichiometry == ((0, 2),)
    assert candidate.score_components["virtual_copy_penalty"] == -0.05
    assert included.self_reactions_considered is True
    assert included.diagnostics.self_reaction_pathway_count == 2
    assert (
        "SELF_REACTION_PATHWAYS_ASSUME_MULTIPLE_EQUIVALENTS_OF_ONE_INPUT"
        in included.warnings
    )
    assert excluded.self_reactions_considered is False
    assert excluded.candidates == ()


def test_transition_metal_profile_favors_structural_cross_coupling_evidence() -> None:
    library = _library(
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    )

    baseline = predict_products("Brc1ccc(N)cc1", library)
    palladium = predict_products(
        "Brc1ccc(N)cc1",
        library,
        condition_profile={
            "strategy": "transition_metal_catalysis",
            "catalyst_family": "palladium",
        },
    )

    evidence = palladium.candidates[0].condition_profile_evidence
    assert palladium.condition_profile_supplied is True
    assert palladium.condition_profile.catalyst_family == "palladium"
    assert evidence.evaluated is True
    assert evidence.score_adjustment == 0.12
    assert evidence.matched_rules == ("transition_metal_cross_coupling",)
    assert palladium.candidates[0].score == round(
        baseline.candidates[0].score + 0.12,
        8,
    )


def test_uncatalyzed_profile_rejects_unactivated_bromoaniline_self_coupling() -> None:
    library = _library(
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    )

    result = predict_products(
        "Brc1ccc(N)cc1",
        library,
        condition_profile={"strategy": "thermal"},
    )

    assert result.status == "no_supported_product"
    assert result.candidates == ()
    assert result.diagnostics.condition_profile_conflict_count == 2
    assert "CONDITION_PROFILE_EXCLUDED_PATHWAYS" in result.warnings

    audit = assess_proposed_step(
        "Brc1ccc(N)cc1",
        "Nc1ccc(Nc2ccc(Br)cc2)cc1",
        library,
        operator_hint=library.operators[0].operator_id,
        condition_profile={"strategy": "thermal"},
    )
    checks = {check.check_id: check.status for check in audit.checks}
    assert audit.targeted_replay_status == "structurally_reproduced"
    assert audit.blind_prediction.candidates == ()
    assert audit.disposition == "condition_incompatible"
    assert audit.validity == "contradicted"
    assert checks["condition_compatibility"] == "warning"


def test_uncatalyzed_profile_retains_activated_aryl_substitution() -> None:
    library = _library(
        "Brc1ccc([N+](=O)[O-])cc1.N>>Nc1ccc([N+](=O)[O-])cc1"
    )

    result = predict_products(
        "Brc1ccc([N+](=O)[O-])cc1.N",
        library,
        condition_profile={"strategy": "thermal"},
    )

    assert tuple(candidate.product_smiles for candidate in result.candidates) == (
        "Nc1ccc([N+](=O)[O-])cc1",
    )
    evidence = result.candidates[0].condition_profile_evidence
    assert evidence.compatible is True
    assert evidence.hard_conflicts == ()
    assert evidence.matched_rules == (
        "thermal_activated_aryl_substitution_allowed",
    )
