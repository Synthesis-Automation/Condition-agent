"""Small orchestration layer for condition recommendation."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Optional, Protocol

from condition_recommender import (
    ChemistRankingPreferences,
    GenericConditionRecommender,
    GenericRecommendationResult,
    build_completion_selection,
    propose_reaction_completion,
)
from condition_recommender.reaction_completion import validate_completion_selections
from reactive_taxonomy import RxnMapperProvider

from .contracts import ConditionRequest, ConditionResponse, ConditionReview
from .rendering import render_recommendation


DEFAULT_INDEX_PATH = Path("datasets/literature/compact/generic_index.sqlite")
DEFAULT_INDEX_CANDIDATES = (
    Path("datasets/literature/full/generic_index.sqlite"),
    DEFAULT_INDEX_PATH,
)


class _Recommender(Protocol):
    def recommend(
        self, reaction_smiles: str, **kwargs: Any
    ) -> GenericRecommendationResult:
        """Return a typed generic recommendation result."""


class _ConditionReviewer(Protocol):
    def review(self, result, settings) -> ConditionReview:
        """Return a bounded post-recommendation review."""


@dataclass(frozen=True)
class ConditionCoworker:
    """Coordinate completion, ranking preferences, and one canonical recommender."""

    recommender: _Recommender
    reviewer: Optional[_ConditionReviewer] = None
    startup_warnings: tuple[str, ...] = ()

    @classmethod
    def from_default(
        cls,
        *,
        use_rxnmapper: bool = False,
        include_review: bool = False,
        reviewer: Optional[_ConditionReviewer] = None,
    ) -> "ConditionCoworker":
        """Select the largest compatible default index without hiding staleness."""

        failures = []
        for candidate in DEFAULT_INDEX_CANDIDATES:
            if not candidate.is_file():
                failures.append(f"{candidate}: not found")
                continue
            try:
                kwargs: dict[str, Any] = {
                    "use_rxnmapper": use_rxnmapper,
                    "include_review": include_review,
                }
                if reviewer is not None:
                    kwargs["reviewer"] = reviewer
                coworker = cls.from_path(candidate, **kwargs)
            except (OSError, ValueError) as exc:
                failures.append(f"{candidate}: {exc}")
                continue
            warnings = ()
            if failures:
                warnings = (
                    "Default index fallback selected "
                    f"{candidate}; skipped " + "; ".join(failures),
                )
            return replace(coworker, startup_warnings=warnings)
        detail = "; ".join(failures) or "no default candidates configured"
        raise ValueError(f"No compatible default recommendation index: {detail}")

    @classmethod
    def from_path(
        cls,
        path: str | Path = DEFAULT_INDEX_PATH,
        *,
        use_rxnmapper: bool = False,
        include_review: bool = False,
        reviewer: Optional[_ConditionReviewer] = None,
    ) -> "ConditionCoworker":
        """Load one validated recommendation index for repeated requests."""

        mapping_provider = None
        if use_rxnmapper:
            if not RxnMapperProvider.is_available():
                raise RuntimeError(
                    "RXNMapper is unavailable; install requirements-mapping.txt"
                )
            mapping_provider = RxnMapperProvider()
        return cls(
            recommender=GenericConditionRecommender.from_path(
                path,
                mapping_provider=mapping_provider,
                include_review=include_review,
            ),
            reviewer=reviewer,
        )

    @staticmethod
    def prepare_reaction(reaction_smiles: str) -> dict[str, Any]:
        """Return a typed source-completion proposal before recommendation."""

        return propose_reaction_completion(reaction_smiles.strip()).to_dict()

    def recommend(self, request: ConditionRequest) -> ConditionResponse:
        """Execute the condition-first flow without LLM chemistry routing."""

        reaction_smiles = request.reaction_smiles.strip()
        selections = ()
        if request.completion_choices:
            proposal = propose_reaction_completion(reaction_smiles)
            expected = {item.requirement_id for item in proposal.requirements}
            supplied = [item.requirement_id for item in request.completion_choices]
            if set(supplied) != expected or len(supplied) != len(expected):
                raise ValueError("REACTION_COMPLETION_CHOICES_INCOMPLETE")
            selections = tuple(
                build_completion_selection(
                    proposal,
                    choice.requirement_id,
                    option_id=choice.option_id,
                    custom_identifier=choice.custom_identifier,
                )
                for choice in request.completion_choices
            )
            validate_completion_selections(proposal, selections)

        preferences = ChemistRankingPreferences(
            profile_id=request.ranking_profile,
            weights=dict(request.ranking_weights),
            customized=bool(request.ranking_weights),
        )
        candidate_pool_size = request.top_k
        if request.review.mode != "off":
            candidate_pool_size = min(
                50,
                max(
                    request.top_k * 2,
                    request.review.max_candidates,
                ),
            )
        result = self.recommender.recommend(
            reaction_smiles,
            top_k=candidate_pool_size,
            minimum_pool_size=request.minimum_pool_size,
            unrestricted_fallback=request.unrestricted_fallback,
            ranking_preferences=preferences,
            completion_selections=selections,
            preferred_reaction_ids=request.preferred_reaction_ids,
        )
        review = None
        if request.review.mode != "off":
            if self.reviewer is None:
                review = ConditionReview(
                    status="failed",
                    provider=request.review.provider,
                    model=request.review.model,
                    warning="LLM review is not configured for this coworker instance",
                    presentation_recipe_ids=tuple(
                        item.recipe_id for item in result.recommendations
                    ),
                )
            else:
                review = self.reviewer.review(result, request.review)
        return ConditionResponse(
            request=request,
            result=result,
            answer=render_recommendation(
                result,
                review,
                display_limit=request.top_k,
            ),
            review=review,
        )
