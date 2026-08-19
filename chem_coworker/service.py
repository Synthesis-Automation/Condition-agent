"""Small orchestration layer for condition recommendation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

from condition_recommender import (
    ChemistRankingPreferences,
    GenericConditionRecommender,
    GenericRecommendationResult,
    build_completion_selection,
    propose_reaction_completion,
)
from condition_recommender.reaction_completion import validate_completion_selections
from reactive_taxonomy import RxnMapperProvider

from .contracts import ConditionRequest, ConditionResponse
from .rendering import render_recommendation


DEFAULT_INDEX_PATH = Path("datasets/literature/full/generic_index.sqlite")


class _Recommender(Protocol):
    def recommend(self, reaction_smiles: str, **kwargs: Any) -> GenericRecommendationResult:
        """Return a typed generic recommendation result."""


@dataclass(frozen=True)
class ConditionCoworker:
    """Coordinate completion, ranking preferences, and one canonical recommender."""

    recommender: _Recommender

    @classmethod
    def from_path(
        cls,
        path: str | Path = DEFAULT_INDEX_PATH,
        *,
        use_rxnmapper: bool = False,
        include_review: bool = False,
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
            GenericConditionRecommender.from_path(
                path,
                mapping_provider=mapping_provider,
                include_review=include_review,
            )
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
        result = self.recommender.recommend(
            reaction_smiles,
            top_k=request.top_k,
            minimum_pool_size=request.minimum_pool_size,
            unrestricted_fallback=request.unrestricted_fallback,
            ranking_preferences=preferences,
            completion_selections=selections,
            preferred_reaction_ids=request.preferred_reaction_ids,
        )
        return ConditionResponse(
            request=request,
            result=result,
            answer=render_recommendation(result),
        )
