"""Deterministic recommendation filtering by confirmed condition constraints."""

from __future__ import annotations

from condition_registry import ConditionConstraintSet, normalize_condition_constraint
from condition_recommender.constraints import apply_condition_constraints
from condition_recommender.models import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)


def _recommendation(rank: int, substance_id: str) -> GenericConditionRecommendation:
    trace = RecommendationScoreTrace(
        similarity_components={},
        similarity_contributions={},
        ranking_components={},
        ranking_contributions={},
        applied_ranking_weights={},
        independent_evidence_count=1,
        observed_outcome_count=0,
        pool_yield_prior_pct=None,
        definition_versions={},
    )
    return GenericConditionRecommendation(
        rank=rank,
        recipe_id=f"recipe:{rank}",
        recipe_core_id=f"core:{rank}",
        recipe_variant_ids=(),
        resolved_recipe={
            "catalysts": [
                {"substance_id": substance_id, "primary_role": "metal_catalyst"}
            ]
        },
        score=1.0 / rank,
        similarity_score=1.0 / rank,
        compatibility_score=1.0,
        expected_yield_pct=None,
        support=1,
        observation_support=1,
        reference_support=1,
        condition_series_support=1,
        dataset_support=1,
        retrieval_level="generic_signature",
        precedent_reaction_ids=(),
        precedent_reference_ids=(),
        explanation=(),
        score_trace=trace,
    )


def test_confirmed_constraint_filters_and_reranks_inside_recommender_layer() -> None:
    constraint = normalize_condition_constraint(
        "excluded_substance",
        "Pd(PPh3)4",
        provenance="confirmed_user",
    ).constraint
    assert constraint is not None
    result = GenericRecommendationResult(
        query_reaction_smiles="CCBr.O>>CCO",
        valid=True,
        recommendations=(
            _recommendation(1, "cas:14221-01-3"),
            _recommendation(2, "cas:7440-50-8"),
        ),
    )

    constrained = apply_condition_constraints(
        result,
        ConditionConstraintSet((constraint,)),
        top_k=5,
    )

    assert [item.recipe_id for item in constrained.recommendations] == ["recipe:2"]
    assert constrained.recommendations[0].rank == 1
    assert constrained.condition_constraint_trace[0]["recipe_id"] == "recipe:1"
    assert result.recommendations[0].rank == 1


def test_all_filtered_result_abstains_explicitly() -> None:
    constraint = normalize_condition_constraint(
        "excluded_substance",
        "Pd(PPh3)4",
        provenance="explicit_user",
    ).constraint
    assert constraint is not None
    result = GenericRecommendationResult(
        query_reaction_smiles="CCBr.O>>CCO",
        valid=True,
        recommendations=(_recommendation(1, "cas:14221-01-3"),),
    )

    constrained = apply_condition_constraints(
        result,
        ConditionConstraintSet((constraint,)),
        top_k=5,
    )

    assert not constrained.valid
    assert constrained.error == "NO_RECOMMENDATIONS_MATCH_CONFIRMED_CONSTRAINTS"
