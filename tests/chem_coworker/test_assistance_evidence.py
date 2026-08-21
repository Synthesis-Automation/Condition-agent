"""Losslessness and immutability tests for condition evidence projection."""

from __future__ import annotations

from condition_recommender.models import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
    RetrievalLevelTrace,
)

from chem_coworker.assistance import project_condition_result


def _result() -> GenericRecommendationResult:
    trace = RecommendationScoreTrace(
        similarity_components={"edit": 1.0, "environment": 0.8},
        similarity_contributions={"edit": 0.5, "environment": 0.3},
        ranking_components={"similarity": 0.8, "yield": 74.0},
        ranking_contributions={"similarity": 0.64, "yield": 0.1},
        applied_ranking_weights={"similarity": 0.8, "yield": 0.2},
        independent_evidence_count=2,
        observed_outcome_count=1,
        pool_yield_prior_pct=65.0,
        definition_versions={"ranking": "ranking.v1@1.0"},
    )
    recommendation = GenericConditionRecommendation(
        rank=1,
        recipe_id="recipe:internal-long-id",
        recipe_core_id="recipe-core:1",
        recipe_variant_ids=("recipe-variant:1",),
        resolved_recipe={
            "components": [
                {"substance_id": "substance:palladium", "role": "catalyst"}
            ]
        },
        score=0.87,
        similarity_score=0.8,
        compatibility_score=1.0,
        expected_yield_pct=74.0,
        support=3,
        observation_support=3,
        reference_support=2,
        condition_series_support=2,
        dataset_support=1,
        retrieval_level="generic_signature",
        precedent_reaction_ids=("reaction:1", "reaction:2"),
        precedent_reference_ids=("reference:1", "reference:2"),
        explanation=("Matched the observed bond edit.",),
        score_trace=trace,
        synthesis_protocol={"temperature_c": 80.0},
        factor_evidence={"edit": {"matched": True}},
        precedent_reaction_smiles=("CCBr.B(O)O>>CCC",),
        precedent_reaction_contexts=({"source": "fixture"},),
        compatibility_evidence=("compatibility.rule.allowed",),
        cautions=("Limited outcome count.",),
    )
    return GenericRecommendationResult(
        query_reaction_smiles="CCBr.B(O)O>>CCC",
        effective_query_reaction_smiles="CCBr.B(O)O>>CCC",
        valid=True,
        query_signature_id="signature:1",
        query_reaction_core_id="core:1",
        query_edit_hypothesis_ids=("edit-hypothesis:1",),
        external_mapping_status="accepted",
        external_mapping_provider="fixture-mapper",
        external_mapping_confidence=0.95,
        reaction_label={"label": "optional"},
        named_family="suzuki",
        transformation_class="coupling",
        spectator_groups=({"component": "base"},),
        reaction_partners=({"role": "electrophile"}, {"role": "partner"}),
        ranking_preferences={"profile_id": "default"},
        retrieval_definition_version="retrieval.v1@1.0",
        retrieval_level="generic_signature",
        candidate_count=5,
        independent_candidate_count=4,
        compatible_candidate_count=3,
        independent_compatible_candidate_count=2,
        excluded_candidate_count=2,
        retrieval_trace=(
            RetrievalLevelTrace(
                level="generic_signature",
                candidate_count=5,
                independent_candidate_count=4,
                compatible_candidate_count=3,
                independent_compatible_candidate_count=2,
                excluded_candidate_count=2,
                minimum_independent_support=1,
                status="selected",
            ),
        ),
        recommendations=(recommendation,),
        warnings=("Mapping and reconstruction agree.",),
    )


def test_projection_is_read_only_and_uses_request_local_aliases() -> None:
    result = _result()
    before = result.to_dict()

    projection = project_condition_result(result)

    assert result.to_dict() == before
    assert projection.resolve_candidate("candidate-1") == "recipe:internal-long-id"
    assert projection.initial_packet()["candidate_aliases"] == ["candidate-1"]


def test_candidate_inspection_preserves_full_score_and_support_evidence() -> None:
    projection = project_condition_result(_result())
    packet = projection.inspection_packet("candidate-1")
    evidence = {item["evidence_id"]: item["payload"] for item in packet["evidence"]}

    ranking = evidence["candidate-1.ranking"]
    support = evidence["candidate-1.support"]
    recipe = evidence["candidate-1.recipe"]
    assert ranking["score_trace"]["similarity_components"]["edit"] == 1.0
    assert ranking["factor_evidence"]["edit"]["matched"] is True
    assert support["observation_support"] == 3
    assert support["precedent_reaction_smiles"] == ["CCBr.B(O)O>>CCC"]
    assert recipe["resolved_recipe"]["components"][0]["role"] == "catalyst"


def test_projection_includes_mapping_retrieval_and_uncertainty() -> None:
    projection = project_condition_result(_result())
    packet = projection.initial_packet()
    evidence = {item["evidence_id"]: item for item in packet["evidence"]}

    assert evidence["query.mapping"]["payload"]["confidence"] == 0.95
    assert evidence["query.retrieval"]["payload"]["retrieval_trace"][0][
        "status"
    ] == "selected"
    assert evidence["candidate-1.summary"]["uncertainty"] == (
        "Limited outcome count."
    )


def test_unknown_aliases_and_single_candidate_comparisons_fail_closed() -> None:
    projection = project_condition_result(_result())

    try:
        projection.inspection_packet("recipe:internal-long-id")
    except ValueError as exc:
        assert "unknown condition candidate alias" in str(exc)
    else:
        raise AssertionError("internal IDs must not be accepted as aliases")

    try:
        projection.comparison_packet(("candidate-1",))
    except ValueError as exc:
        assert "at least two" in str(exc)
    else:
        raise AssertionError("single candidate comparison must fail")

