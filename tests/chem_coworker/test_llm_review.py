"""Tests for the bounded post-recommendation LLM review layer."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from condition_recommender import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)

from chem_coworker.contracts import ConditionReviewSettings
from chem_coworker.review import (
    CandidateReviewPayload,
    ConditionGroupPayload,
    ConditionReviewPayload,
    LLMConditionReviewer,
    ReviewTransportResult,
)


def _recommendation(rank: int, score: float, support: int = 2):
    return GenericConditionRecommendation(
        rank=rank,
        recipe_id=f"recipe-{rank}",
        recipe_core_id=f"core-{rank}",
        recipe_variant_ids=(),
        resolved_recipe={"components": [{"canonical_name": f"base {rank}"}]},
        score=score,
        similarity_score=score,
        compatibility_score=1.0,
        expected_yield_pct=None,
        support=support,
        observation_support=support,
        reference_support=support,
        condition_series_support=support,
        dataset_support=1,
        retrieval_level="reaction_facet_exact",
        precedent_reaction_ids=(f"rxn-{rank}",),
        precedent_reference_ids=(f"ref-{rank}",),
        explanation=("matching structural reaction facet",),
        score_trace=RecommendationScoreTrace(
            similarity_components={},
            similarity_contributions={},
            ranking_components={},
            ranking_contributions={},
            applied_ranking_weights={},
            independent_evidence_count=support,
            observed_outcome_count=support,
            pool_yield_prior_pct=None,
            definition_versions={},
        ),
    )


@dataclass
class FakeTransport:
    payload: ConditionReviewPayload
    packet: Mapping[str, Any] | None = None

    def complete(self, evidence_packet, settings):
        self.packet = evidence_packet
        return ReviewTransportResult(
            payload=self.payload,
            response_id="resp-test",
            input_tokens=100,
            output_tokens=20,
        )


def test_auto_review_skips_single_high_confidence_recommendation() -> None:
    result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
        retrieval_level="reaction_facet_exact",
        recommendations=(_recommendation(1, 0.9),),
    )
    reviewer = LLMConditionReviewer()

    review = reviewer.review(
        result,
        ConditionReviewSettings(mode="auto"),
    )

    assert review.status == "skipped"
    assert review.trigger_reasons == ("skipped_no_uncertainty_signal",)


def test_review_can_reorder_presentation_without_mutating_domain_ranking() -> None:
    first = _recommendation(1, 0.91)
    second = _recommendation(2, 0.89)
    payload = ConditionReviewPayload(
        summary="The first recipe has a selectivity concern.",
        candidates=[
            CandidateReviewPayload(
                recipe_id="recipe-1",
                verdict="flag",
                issue_codes=["chemoselectivity"],
                evidence_ids=["candidate.recipe-1"],
                rationale="Two plausible reactive sites are not distinguished.",
                confidence=0.8,
            ),
            CandidateReviewPayload(
                recipe_id="recipe-2",
                verdict="keep",
                issue_codes=[],
                evidence_ids=["candidate.recipe-2"],
                rationale="No material conflict is visible in the supplied evidence.",
                confidence=0.7,
            ),
        ],
        groups=[
            ConditionGroupPayload(
                group_id="flagged-strategy",
                member_recipe_ids=["recipe-1"],
                grouping_basis=["same_strategy"],
                evidence_ids=["candidate.recipe-1"],
                rationale="This recipe is a distinct strategy.",
            ),
            ConditionGroupPayload(
                group_id="retained-strategy",
                member_recipe_ids=["recipe-2"],
                grouping_basis=["same_strategy"],
                evidence_ids=["candidate.recipe-2"],
                rationale="This recipe is a distinct strategy.",
            ),
        ],
        questions=["Which amine should react?"],
    )
    transport = FakeTransport(payload)
    reviewer = LLMConditionReviewer(transport)
    result = GenericRecommendationResult(
        query_reaction_smiles="CN.CN.ClC>>product",
        valid=True,
        retrieval_level="reaction_facet_exact",
        recommendations=(first, second),
    )

    review = reviewer.review(
        result,
        ConditionReviewSettings(mode="always"),
    )

    assert review.status == "completed"
    assert review.presentation_recipe_ids == ("recipe-2", "recipe-1")
    assert [item.recipe_id for item in result.recommendations] == [
        "recipe-1",
        "recipe-2",
    ]
    assert transport.packet is not None
    assert transport.packet["review_candidate_ids"] == ("recipe-1", "recipe-2")
    assert review.response_id == "resp-test"


def test_unknown_model_evidence_reference_fails_closed() -> None:
    payload = ConditionReviewPayload(
        summary="Concern found.",
        candidates=[
            CandidateReviewPayload(
                recipe_id="recipe-1",
                verdict="downrank",
                issue_codes=["other"],
                evidence_ids=["invented.evidence"],
                rationale="Unsupported rationale.",
                confidence=0.9,
            )
        ],
        groups=[
            ConditionGroupPayload(
                group_id="strategy-1",
                member_recipe_ids=["recipe-1"],
                grouping_basis=["same_strategy"],
                evidence_ids=["candidate.recipe-1"],
                rationale="This recipe is a distinct strategy.",
            )
        ],
        questions=[],
    )
    result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
        recommendations=(_recommendation(1, 0.9),),
    )

    review = LLMConditionReviewer(FakeTransport(payload)).review(
        result,
        ConditionReviewSettings(mode="always"),
    )

    assert review.status == "failed"
    assert review.presentation_recipe_ids == ("recipe-1",)
    assert "unknown evidence IDs" in (review.warning or "")


def test_similar_recipe_variants_collapse_to_one_strategy_representative() -> None:
    payload = ConditionReviewPayload(
        summary="Two carbonate/aqueous solvent recipes share one catalyst strategy.",
        candidates=[
            CandidateReviewPayload(
                recipe_id=f"recipe-{rank}",
                verdict="keep",
                issue_codes=[],
                evidence_ids=[f"candidate.recipe-{rank}"],
                rationale="The supplied evidence supports this recipe.",
                confidence=0.8,
            )
            for rank in (1, 2, 3)
        ],
        groups=[
            ConditionGroupPayload(
                group_id="pd-pph3-carbonate",
                member_recipe_ids=["recipe-1", "recipe-2"],
                grouping_basis=[
                    "catalyst_system",
                    "base_family",
                    "solvent_system",
                ],
                evidence_ids=["candidate.recipe-1", "candidate.recipe-2"],
                rationale=(
                    "The same palladium catalyst strategy differs only in carbonate "
                    "base and aqueous ethereal solvent."
                ),
            ),
            ConditionGroupPayload(
                group_id="pd-dppf",
                member_recipe_ids=["recipe-3"],
                grouping_basis=["catalyst_system"],
                evidence_ids=["candidate.recipe-3"],
                rationale="The catalyst/ligand system is materially different.",
            ),
        ],
        questions=[],
    )
    result = GenericRecommendationResult(
        query_reaction_smiles="B(O)O.ClC>>CC",
        valid=True,
        recommendations=tuple(
            _recommendation(rank, score)
            for rank, score in ((1, 0.72), (2, 0.71), (3, 0.69))
        ),
    )

    review = LLMConditionReviewer(FakeTransport(payload)).review(
        result,
        ConditionReviewSettings(mode="always"),
    )

    assert review.status == "completed"
    assert review.presentation_recipe_ids == ("recipe-1", "recipe-3")
    assert review.groups[0].representative_recipe_id == "recipe-1"
    assert review.groups[0].member_recipe_ids == ("recipe-1", "recipe-2")
