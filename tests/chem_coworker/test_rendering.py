"""Rendering tests for evidence-bearing recipe summaries."""

from condition_recommender import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)

from chem_coworker.rendering import render_recommendation


def test_renderer_includes_recipe_evidence_cautions_and_precedents() -> None:
    score_trace = RecommendationScoreTrace(
        similarity_components={},
        similarity_contributions={},
        ranking_components={},
        ranking_contributions={},
        applied_ranking_weights={},
        independent_evidence_count=2,
        observed_outcome_count=2,
        pool_yield_prior_pct=70.0,
        definition_versions={},
    )
    recommendation = GenericConditionRecommendation(
        rank=1,
        recipe_id="recipe-1",
        recipe_core_id="core-1",
        recipe_variant_ids=("variant-1",),
        resolved_recipe={
            "components": [
                {"canonical_name": "palladium acetate", "primary_role": "catalyst"},
                {"canonical_name": "potassium carbonate", "primary_role": "base"},
            ]
        },
        score=0.91,
        similarity_score=0.88,
        compatibility_score=1.0,
        expected_yield_pct=74.5,
        support=2,
        observation_support=2,
        reference_support=2,
        condition_series_support=2,
        dataset_support=1,
        retrieval_level="bond_edit_exact",
        precedent_reaction_ids=("rxn-1", "rxn-2"),
        precedent_reference_ids=("ref-1", "ref-2"),
        explanation=("matching C-N bond edit",),
        score_trace=score_trace,
        cautions=("substrate is more hindered",),
    )
    result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
        recommendations=(recommendation,),
        retrieval_level="bond_edit_exact",
        candidate_count=2,
    )

    text = render_recommendation(result)

    assert "palladium acetate (catalyst)" in text
    assert "matching C-N bond edit" in text
    assert "substrate is more hindered" in text
    assert "rxn-1, rxn-2" in text
