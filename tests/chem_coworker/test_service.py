"""Behavior tests for the condition-first orchestration service."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from condition_recommender import GenericRecommendationResult

from chem_coworker import (
    ConditionCoworker,
    ConditionRequest,
    ConditionReview,
    ConditionReviewSettings,
)
import chem_coworker.service as service_module


@dataclass
class FakeRecommender:
    result: GenericRecommendationResult
    calls: list[tuple[str, dict[str, Any]]] = field(default_factory=list)

    def recommend(
        self, reaction_smiles: str, **kwargs: Any
    ) -> GenericRecommendationResult:
        self.calls.append((reaction_smiles, kwargs))
        return self.result


@dataclass
class FakeReviewer:
    calls: list[tuple[GenericRecommendationResult, ConditionReviewSettings]] = field(
        default_factory=list
    )

    def review(self, result, settings) -> ConditionReview:
        self.calls.append((result, settings))
        return ConditionReview(
            status="completed",
            provider=settings.provider,
            model=settings.model,
            summary="No contextual conflict found.",
        )


def test_recommendation_is_delegated_once_to_generic_system() -> None:
    domain_result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
        named_family=None,
        transformation_class="generic_c_n_coupling",
        retrieval_level="bond_edit_exact",
        candidate_count=3,
        warnings=("LOW_INDEPENDENT_SUPPORT",),
    )
    recommender = FakeRecommender(domain_result)
    coworker = ConditionCoworker(recommender)
    request = ConditionRequest(
        reaction_smiles="  C.N>>CN  ",
        top_k=4,
        ranking_profile="yield_focused",
        preferred_reaction_ids=("rxn-7",),
    )

    response = coworker.recommend(request)

    assert response.result is domain_result
    assert len(recommender.calls) == 1
    reaction, arguments = recommender.calls[0]
    assert reaction == "C.N>>CN"
    assert arguments["top_k"] == 4
    assert arguments["ranking_preferences"].profile_id == "yield_focused"
    assert arguments["preferred_reaction_ids"] == ("rxn-7",)
    assert "family: unassigned (generic structural retrieval)" in response.answer
    assert "LOW_INDEPENDENT_SUPPORT" in response.answer


def test_invalid_result_abstains_and_preserves_full_warning() -> None:
    long_warning = "EVIDENCE:" + "x" * 4_000
    result = GenericRecommendationResult(
        query_reaction_smiles="invalid",
        valid=False,
        error="INVALID_REACTION",
        warnings=(long_warning,),
    )
    coworker = ConditionCoworker(FakeRecommender(result))

    response = coworker.recommend(ConditionRequest(reaction_smiles="invalid"))
    payload = response.to_dict()

    assert not response.valid
    assert response.answer.startswith("No condition recommendation")
    assert "Recommended conditions:" not in response.answer
    assert payload["result"]["warnings"] == (long_warning,)
    assert len(payload["result"]["warnings"][0]) > 3_000


def test_request_rejects_unsafe_ranking_weights() -> None:
    try:
        ConditionRequest(
            reaction_smiles="C.N>>CN",
            ranking_weights={"expected_yield": -1.0},
        )
    except ValueError as exc:
        assert str(exc) == "ranking weights must be non-negative"
    else:
        raise AssertionError("negative ranking weights were accepted")


def test_optional_llm_review_is_applied_after_deterministic_recommendation() -> None:
    domain_result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
    )
    recommender = FakeRecommender(domain_result)
    reviewer = FakeReviewer()
    coworker = ConditionCoworker(recommender, reviewer=reviewer)

    response = coworker.recommend(
        ConditionRequest(
            reaction_smiles="C.N>>CN",
            review=ConditionReviewSettings(mode="always"),
        )
    )

    assert response.result is domain_result
    assert response.review is not None
    assert response.review.status == "completed"
    assert reviewer.calls == [(domain_result, response.request.review)]
    assert recommender.calls[0][1]["top_k"] == 10
    assert "No contextual conflict found" in response.answer


def test_default_index_selection_skips_stale_full_artifact(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stale = tmp_path / "full.sqlite"
    current = tmp_path / "compact.sqlite"
    stale.touch()
    current.touch()
    domain_result = GenericRecommendationResult(
        query_reaction_smiles="C.N>>CN",
        valid=True,
    )

    def fake_from_path(
        cls,
        path,
        *,
        use_rxnmapper=False,
        include_review=False,
    ):
        del use_rxnmapper, include_review
        if Path(path) == stale:
            raise ValueError("Unsupported generic reaction index schema")
        return cls(FakeRecommender(domain_result))

    monkeypatch.setattr(
        service_module,
        "DEFAULT_INDEX_CANDIDATES",
        (stale, current),
    )
    monkeypatch.setattr(
        ConditionCoworker,
        "from_path",
        classmethod(fake_from_path),
    )

    coworker = ConditionCoworker.from_default()

    assert len(coworker.startup_warnings) == 1
    assert str(current) in coworker.startup_warnings[0]
    assert "Unsupported generic reaction index schema" in coworker.startup_warnings[0]
