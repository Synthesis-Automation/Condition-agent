"""Application orchestration for chemistry-first one-step retrosynthesis."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Protocol, Sequence

from condition_recommender.processed_release import resolve_processed_release

from core_retrosynthesis import (
    GenericTemplateLibrary,
    StrategyProposal,
    disconnect_operator_ladder,
    disconnect_operator_ladder_detailed,
    group_strategy_candidates,
    load_generic_library,
    load_retrosynthesis_ranking_policy,
    recommend_retrosynthesis_conditions,
)

from .contracts import (
    RetrosynthesisRequest,
    RetrosynthesisResponse,
    RetrosynthesisReview,
    RetrosynthesisStrategyCondition,
)
from .retrosynthesis_rendering import render_retrosynthesis
from reactive_taxonomy.disconnection_focus import prepare_disconnection_focus


_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RETROSYNTHESIS_LIBRARY_PATH = _REPOSITORY_ROOT / Path(
    "datasets/processed_datasets/manifest.json"
)


class _RetrosynthesisReviewer(Protocol):
    def review(
        self,
        strategies: Sequence[StrategyProposal],
        condition_evidence: Sequence[RetrosynthesisStrategyCondition],
        settings,
    ) -> RetrosynthesisReview:
        """Return bounded annotations for existing strategies."""


class _ConditionRecommender(Protocol):
    def recommend(self, reaction_smiles: str, **kwargs):
        """Return deterministic condition evidence."""


@dataclass(frozen=True)
class RetrosynthesisCoworker:
    """Coordinate validated strategy search, condition evidence, and review."""

    library: GenericTemplateLibrary
    library_path: Optional[Path] = None
    condition_recommender: Optional[_ConditionRecommender] = None
    reviewer: Optional[_RetrosynthesisReviewer] = None
    startup_warnings: tuple[str, ...] = ()

    @classmethod
    def from_default(
        cls,
        *,
        condition_recommender: Optional[_ConditionRecommender] = None,
        reviewer: Optional[_RetrosynthesisReviewer] = None,
    ) -> "RetrosynthesisCoworker":
        """Pin the published complete corpus's validated operator library."""
        return cls.from_path(
            DEFAULT_RETROSYNTHESIS_LIBRARY_PATH,
            condition_recommender=condition_recommender,
            reviewer=reviewer,
        )

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        condition_recommender: Optional[_ConditionRecommender] = None,
        reviewer: Optional[_RetrosynthesisReviewer] = None,
    ) -> "RetrosynthesisCoworker":
        """Load one validated generic operator library."""

        resolved = Path(path).expanduser().resolve()
        if resolved.is_dir() or resolved.name == "manifest.json":
            resolved = resolve_processed_release(resolved).artifact("retro_library")
        return cls(
            library=load_generic_library(resolved),
            library_path=resolved,
            condition_recommender=condition_recommender,
            reviewer=reviewer,
        )

    def disconnect(self, request: RetrosynthesisRequest) -> RetrosynthesisResponse:
        """Return graph-validated strategies with bounded advisory review."""

        target = request.target_smiles.strip()
        pool_size = request.top_k
        if request.review.mode != "off":
            pool_size = min(
                50,
                max(request.top_k * 2, request.review.max_candidates),
            )
        try:
            focus = None
            search_diagnostics = None
            if request.required_disconnection_bond is not None:
                focus = prepare_disconnection_focus(
                    target, request.required_disconnection_bond, request.focus_target_smiles,
                )
            elif request.focus_target_smiles is not None:
                raise ValueError("focus_target_smiles requires required_disconnection_bond")
            policy = load_retrosynthesis_ranking_policy()
            flat_pool_size = min(
                request.max_candidates_to_validate,
                max(
                    pool_size,
                    pool_size * policy.candidate_pool_multiplier,
                    pool_size * request.max_realizations_per_strategy,
                ),
            )
            options = dict(top_k=flat_pool_size,
                           max_templates_to_apply=request.max_templates_to_apply,
                           max_candidates_to_validate=request.max_candidates_to_validate,
                           use_context=request.use_context, include_l0=request.include_l0)
            if focus is not None:
                candidates, diagnostics = disconnect_operator_ladder_detailed(
                    target, self.library, **options,
                    required_disconnection_bond=focus.atom_ids,
                    focus_target_smiles=focus.canonical_target_smiles,
                )
                search_diagnostics = diagnostics.to_dict()
            else:
                candidates = disconnect_operator_ladder(target, self.library, **options)
            identity_candidates = tuple(
                candidate for candidate in candidates if candidate.strategy_id
            )
            strategies = group_strategy_candidates(
                identity_candidates,
                top_k_strategies=pool_size,
                max_realizations_per_strategy=request.max_realizations_per_strategy,
            )
        except ValueError as exc:
            response = RetrosynthesisResponse(
                request=request,
                valid=False,
                error=str(exc),
                library_path=str(self.library_path) if self.library_path else None,
            )
            return replace(response, answer=render_retrosynthesis(response))

        condition_evidence = self._condition_evidence(strategies, request)
        review = self._review(strategies, condition_evidence, request)
        warnings = []
        if focus is not None:
            warnings.append(
                "Focused search retains ordinary template retrieval and per-level budgets; "
                "no qualifying result is a bounded search outcome, not synthetic impossibility."
            )
        omitted_identity_count = len(candidates) - len(identity_candidates)
        if omitted_identity_count:
            warnings.append(
                f"{omitted_identity_count} validated transformations without a "
                "complete target-site strategy identity were omitted from the "
                "one-step disconnection view."
            )
        if not strategies:
            warnings.append(
                "No forward-validated single-step strategy was found in the "
                "loaded operator library."
            )
        if request.include_conditions and self.condition_recommender is None:
            warnings.append(
                "Condition evidence was requested but no condition recommender "
                "is configured."
            )
        response = RetrosynthesisResponse(
            request=request,
            valid=True,
            strategies=tuple(strategies),
            condition_evidence=condition_evidence,
            review=review,
            bond_focus=focus,
            search_diagnostics=search_diagnostics,
            warnings=tuple(warnings),
            library_path=str(self.library_path) if self.library_path else None,
        )
        return replace(response, answer=render_retrosynthesis(response))

    def _condition_evidence(
        self,
        strategies: Sequence[StrategyProposal],
        request: RetrosynthesisRequest,
    ) -> tuple[RetrosynthesisStrategyCondition, ...]:
        if not request.include_conditions or self.condition_recommender is None:
            return ()
        values = []
        for strategy in strategies:
            value = self.recommend_strategy_conditions(strategy, request)
            if value is not None:
                values.append(value)
        return tuple(values)

    def recommend_strategy_conditions(
        self,
        strategy: StrategyProposal,
        request: RetrosynthesisRequest,
    ) -> RetrosynthesisStrategyCondition | None:
        """Recommend conditions for one existing validated representative."""

        if not request.include_conditions or self.condition_recommender is None:
            return None
        candidate = strategy.representative
        reaction_smiles = (
            candidate.condition_query_reaction_smiles
            or candidate.proposed_reaction_smiles
        )
        evidence = recommend_retrosynthesis_conditions(
            reaction_smiles,
            self.condition_recommender,
            condition_top_k=request.condition_top_k,
            minimum_pool_size=request.condition_minimum_pool_size,
            unrestricted_fallback=request.unrestricted_condition_fallback,
            preferred_reaction_ids=candidate.precedent_reaction_ids,
        )
        return RetrosynthesisStrategyCondition(
            strategy_id=strategy.strategy_id,
            evidence=evidence,
        )

    def _review(
        self,
        strategies: Sequence[StrategyProposal],
        condition_evidence: Sequence[RetrosynthesisStrategyCondition],
        request: RetrosynthesisRequest,
    ) -> RetrosynthesisReview | None:
        if request.review.mode == "off":
            return None
        if self.reviewer is None:
            return RetrosynthesisReview(
                status="failed",
                provider=request.review.provider,
                model=request.review.model,
                warning="LLM review is not configured for retrosynthesis",
                presentation_strategy_ids=tuple(
                    item.strategy_id for item in strategies
                ),
            )
        return self.reviewer.review(strategies, condition_evidence, request.review)


__all__ = [
    "DEFAULT_RETROSYNTHESIS_LIBRARY_PATH",
    "RetrosynthesisCoworker",
]
