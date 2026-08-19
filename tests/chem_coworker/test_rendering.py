"""Rendering tests for evidence-bearing recipe summaries."""

from condition_recommender import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)

from chem_coworker.contracts import ConditionGroupReview, ConditionReview
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


def test_renderer_outputs_group_representative_and_hides_minor_variants() -> None:
    from dataclasses import replace

    score_trace = RecommendationScoreTrace(
        similarity_components={},
        similarity_contributions={},
        ranking_components={},
        ranking_contributions={},
        applied_ranking_weights={},
        independent_evidence_count=2,
        observed_outcome_count=2,
        pool_yield_prior_pct=None,
        definition_versions={},
    )

    def recommendation(rank: int, name: str) -> GenericConditionRecommendation:
        return GenericConditionRecommendation(
            rank=rank,
            recipe_id=f"recipe-{rank}",
            recipe_core_id=f"core-{rank}",
            recipe_variant_ids=(),
            resolved_recipe={
                "components": [{"canonical_name": name, "primary_role": "catalyst"}]
            },
            score=0.8 - rank / 100,
            similarity_score=0.8,
            compatibility_score=1.0,
            expected_yield_pct=None,
            support=2,
            observation_support=2,
            reference_support=2,
            condition_series_support=2,
            dataset_support=1,
            retrieval_level="reaction_facet_exact",
            precedent_reaction_ids=(),
            precedent_reference_ids=(),
            explanation=(),
            score_trace=score_trace,
        )

    first = recommendation(1, "Pd(PPh3)4 / sodium carbonate")
    second = replace(
        recommendation(2, "Pd(PPh3)4 / potassium carbonate"),
        recipe_id="recipe-2",
    )
    third = recommendation(3, "Pd(dppf)Cl2")
    result = GenericRecommendationResult(
        query_reaction_smiles="B(O)O.ClC>>CC",
        valid=True,
        recommendations=(first, second, third),
    )
    review = ConditionReview(
        status="completed",
        provider="openai",
        model="gpt-5.6-terra",
        summary="Two distinct condition strategies remain.",
        groups=(
            ConditionGroupReview(
                group_id="pd-pph3",
                representative_recipe_id="recipe-1",
                member_recipe_ids=("recipe-1", "recipe-2"),
                grouping_basis=("catalyst_system", "base_family"),
                evidence_ids=("candidate.recipe-1", "candidate.recipe-2"),
                rationale="The catalyst strategy is the same.",
            ),
            ConditionGroupReview(
                group_id="pd-dppf",
                representative_recipe_id="recipe-3",
                member_recipe_ids=("recipe-3",),
                grouping_basis=("catalyst_system",),
                evidence_ids=("candidate.recipe-3",),
                rationale="The ligand system is distinct.",
            ),
        ),
        presentation_recipe_ids=("recipe-1", "recipe-3"),
    )

    text = render_recommendation(result, review, display_limit=5)

    assert "Pd(PPh3)4 / sodium carbonate" in text
    assert "Pd(PPh3)4 / potassium carbonate" not in text
    assert "Pd(dppf)Cl2" in text
    assert "Grouped strategy: 2 recipe variants" in text
