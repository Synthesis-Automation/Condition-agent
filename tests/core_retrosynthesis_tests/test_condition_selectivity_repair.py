"""Precedent-backed condition-selectivity repair regressions."""

from __future__ import annotations

from types import SimpleNamespace

from core_retrosynthesis import (
    assess_condition_selectivity_repairs,
    detect_functional_group_competition,
    load_condition_selectivity_repair_policy,
)


S_ALKYLATION = "Cc1[nH]cnc1CCl.NCCS>>NCCSCc1c(C)[nH]cn1"
N_ALKYLATION = "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1"
NEUTRAL_RECIPE = {
    "medium": "neutral",
    "salt_state": "free_base",
    "temperature_c": 25.0,
}


def _recommendation(
    reactions: tuple[str, ...],
    *,
    compatibility_score: float = 1.0,
) -> dict:
    return {
        "rank": 1,
        "recipe_id": "recipe:neutral",
        "recipe_core_id": "recipe-core:neutral",
        "resolved_recipe": NEUTRAL_RECIPE,
        "compatibility_score": compatibility_score,
        "precedent_reaction_contexts": tuple(
            {
                "reaction_smiles": reaction,
                "reference_id": f"reference:{index}",
            }
            for index, reaction in enumerate(reactions, start=1)
        ),
    }


def _step(
    recommendations: tuple[dict, ...],
    *,
    condition_status: str = "recommended_direct",
):
    warning = detect_functional_group_competition(N_ALKYLATION)
    assert warning is not None
    return SimpleNamespace(
        step_id="step:n-versus-s",
        candidate=SimpleNamespace(
            selectivity_warnings=(warning,),
            condition_query_reaction_smiles=N_ALKYLATION,
            proposed_reaction_smiles=N_ALKYLATION,
        ),
        condition_evidence=SimpleNamespace(
            status=condition_status,
            recommendations=recommendations,
        ),
        condition_selectivity_assessment=None,
    )


def test_condition_selectivity_policy_is_versioned_and_validated() -> None:
    policy = load_condition_selectivity_repair_policy()

    assert policy["definition_id"] == "condition_selectivity_repair.v1"
    assert policy["minimum_independent_references"] == 2
    assert policy["require_direct_condition_retrieval"] is True


def test_direct_independent_precedents_make_intended_endpoint_actionable() -> None:
    assessment = assess_condition_selectivity_repairs(
        _step((_recommendation((N_ALKYLATION, N_ALKYLATION)),))
    )[0]

    assert assessment.status == "supported"
    assert assessment.actionable is True
    assert assessment.desired_probability is not None
    assert assessment.desired_probability > 0.7
    assert assessment.probability_margin is not None
    assert assessment.probability_margin > 0.2
    assert assessment.exact_condition_reference_ids == (
        "reference:1",
        "reference:2",
    )


def test_competing_endpoint_precedents_are_explicitly_non_actionable() -> None:
    assessment = assess_condition_selectivity_repairs(
        _step((_recommendation((S_ALKYLATION, S_ALKYLATION)),))
    )[0]

    assert assessment.status == "competing_outcome_favored"
    assert assessment.actionable is False
    assert assessment.probability_margin is not None
    assert assessment.probability_margin < 0.0
    assert assessment.ranked_outcomes[0]["element"] == "S"


def test_mixed_endpoint_precedents_remain_ambiguous() -> None:
    assessment = assess_condition_selectivity_repairs(
        _step(
            (
                _recommendation(
                    (N_ALKYLATION, N_ALKYLATION, S_ALKYLATION, S_ALKYLATION)
                ),
            )
        )
    )[0]

    assert assessment.status == "ambiguous"
    assert assessment.actionable is False
    assert assessment.normalized_entropy is not None
    assert assessment.normalized_entropy > 0.9


def test_condition_conflict_blocks_otherwise_supported_recipe() -> None:
    assessment = assess_condition_selectivity_repairs(
        _step(
            (
                _recommendation(
                    (N_ALKYLATION, N_ALKYLATION),
                    compatibility_score=0.2,
                ),
            )
        )
    )[0]

    assert assessment.status == "condition_conflict"
    assert assessment.actionable is False


def test_missing_or_fallback_condition_evidence_cannot_trigger_repair() -> None:
    missing = assess_condition_selectivity_repairs(_step(()))[0]
    fallback = assess_condition_selectivity_repairs(
        _step(
            (_recommendation((N_ALKYLATION, N_ALKYLATION)),),
            condition_status="recommended_fallback",
        )
    )[0]

    assert missing.status == "no_condition_evidence"
    assert missing.actionable is False
    assert fallback.status == "insufficient_precedent"
    assert fallback.actionable is False

