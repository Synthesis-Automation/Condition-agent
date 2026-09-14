"""Strategy-grouped single-step retrosynthesis over validated candidates."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Callable, Iterable

from cas_tools import PrecursorRealismAssessment

from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    OperatorLadderDiagnostics,
    StrategyProposal,
)
from .generic_search import (
    _attach_precursor_realism,
    _apply_strategic_candidate_reserve,
    disconnect_generic_target_detailed,
    rank_operator_site_diverse,
    rank_precursor_realism,
)
from .hierarchical_ranking import (
    build_completion_prior_index,
    rank_hierarchical_candidates,
)
from .ranking_policy import load_retrosynthesis_ranking_policy


def group_strategy_candidates(
    candidates: Iterable[GenericDisconnectionCandidate],
    *,
    top_k_strategies: int = 10,
    max_realizations_per_strategy: int = 3,
) -> tuple[StrategyProposal, ...]:
    """Group an ordered validated candidate stream by ``STRAT1`` identity.

    Input order is the existing chemistry-first ranking order.  The first
    realization encountered represents the strategy, so grouping does not
    introduce a second opaque score or change candidate ranking semantics.
    Independent support is aggregated conservatively with ``max`` rather than
    summed because multiple templates may encode overlapping evidence.
    """

    if top_k_strategies < 1:
        raise ValueError("top-k strategies must be positive")
    if max_realizations_per_strategy < 1:
        raise ValueError("maximum realizations per strategy must be positive")

    groups: dict[str, list[GenericDisconnectionCandidate]] = {}
    seen_realizations: dict[str, set[str]] = {}
    for candidate in candidates:
        if candidate.forward_validation_status != "verified_signature":
            raise ValueError("strategy grouping requires verified-signature candidates")
        if not candidate.strategy_id:
            raise ValueError("strategy grouping requires complete STRAT1 identity")
        group = groups.setdefault(candidate.strategy_id, [])
        if group and group[0].target_smiles != candidate.target_smiles:
            raise ValueError("strategy grouping cannot mix target molecules")
        seen = seen_realizations.setdefault(candidate.strategy_id, set())
        realization_key = candidate.precursor_smiles
        if realization_key in seen:
            continue
        seen.add(realization_key)
        group.append(candidate)

    proposals = []
    for strategy_rank, (strategy_id, group) in enumerate(
        tuple(groups.items())[:top_k_strategies],
        start=1,
    ):
        representative = group[0]
        retained = group[:max_realizations_per_strategy]
        proposals.append(
            StrategyProposal(
                strategy_rank=strategy_rank,
                strategy_id=strategy_id,
                operator_id=representative.operator_id,
                disconnection_site_key=representative.disconnection_site_key,
                synthon_signature=representative.synthon_signature,
                representative=representative,
                alternate_realizations=tuple(retained[1:]),
                total_realization_count=len(group),
                independent_reference_support=max(
                    candidate.independent_reference_support for candidate in group
                ),
                precedent_reaction_ids=tuple(
                    sorted(
                        {
                            reaction_id
                            for candidate in group
                            for reaction_id in candidate.precedent_reaction_ids
                            if reaction_id
                        }
                    )
                ),
            )
        )
    return tuple(proposals)


@dataclass(frozen=True)
class StrategySearchResult:
    """Bounded strategy search with explicit effort and shortfall evidence."""

    strategies: tuple[StrategyProposal, ...]
    diagnostics: OperatorLadderDiagnostics
    requested_strategy_count: int
    max_templates_per_level: int
    max_validations_per_level: int
    max_realizations_per_strategy: int
    unresolved_candidates: tuple[GenericDisconnectionCandidate, ...] = ()
    definition_id: str = "single_step_strategy_search.v2"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the canonical grouped result and its search limits."""

        diagnostics = self.diagnostics.to_dict()
        limited = any(
            item.template_budget_excluded_count or item.validation_budget_excluded_count
            for _, item in self.diagnostics.level_diagnostics
        )
        return {
            "definition_id": self.definition_id,
            "strategies": [strategy.to_dict() for strategy in self.strategies],
            "strategy_count": len(self.strategies),
            "requested_strategy_count": self.requested_strategy_count,
            "returned_realization_count": sum(
                len(s.realizations) for s in self.strategies
            ),
            "unresolved_candidates": [
                candidate.to_dict() for candidate in self.unresolved_candidates
            ],
            "search_diagnostics": {
                **diagnostics,
                "budget_limited": limited,
                "strategy_target_met": len(self.strategies)
                >= self.requested_strategy_count,
                "max_templates_per_level": self.max_templates_per_level,
                "max_validations_per_level": self.max_validations_per_level,
                "max_realizations_per_strategy": self.max_realizations_per_strategy,
                "scheduling": "operator_site_round_robin_then_verified_strategy_grouping",
                "incomplete_strategy_identity_count": len(self.unresolved_candidates),
            },
        }


def disconnect_strategies_detailed(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    top_k_strategies: int = 10,
    max_realizations_per_strategy: int = 3,
    max_templates_to_apply: int = 500,
    max_candidates_to_validate: int = 100,
    use_context: bool = True,
    include_l0: bool = True,
    diversify: bool = True,
    use_hierarchical_ranking: bool = True,
    precursor_realism_scorer: (
        Callable[[str], tuple[PrecursorRealismAssessment, ...]] | None
    ) = None,
) -> StrategySearchResult:
    """Search specificity tiers until enough verified strategies are found.

    Operators share template and validation budgets before final grouping.
    Source correspondence and hard graph validation remain authoritative.
    Limits apply per attempted specificity tier, not to the entire request.
    """

    if top_k_strategies < 1:
        raise ValueError("top-k strategies must be positive")
    if max_realizations_per_strategy < 1:
        raise ValueError("maximum realizations per strategy must be positive")
    if max_templates_to_apply < 1:
        raise ValueError("maximum templates to apply must be positive")
    if max_candidates_to_validate < 1:
        raise ValueError("maximum candidates to validate must be positive")

    policy = load_retrosynthesis_ranking_policy()
    prior = (
        build_completion_prior_index(library)
        if diversify and use_hierarchical_ranking
        else None
    )
    candidates = []
    unresolved = {}
    level_diagnostics = []
    strategies: tuple[StrategyProposal, ...] = ()
    for level in ("L2", "L1", "L0") if include_l0 else ("L2", "L1"):
        batch, diagnostics = disconnect_generic_target_detailed(
            target_smiles,
            library,
            levels=(level,),
            top_k=max_candidates_to_validate,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context,
            balance_operator_budget=True,
        )
        level_diagnostics.append((level, diagnostics))
        if precursor_realism_scorer is not None:
            batch = _attach_precursor_realism(batch, precursor_realism_scorer)
        if diversify:
            batch = rank_operator_site_diverse(batch, policy=policy)
            if use_hierarchical_ranking:
                batch = rank_hierarchical_candidates(
                    batch,
                    library,
                    structural_policy=policy,
                    prior_index=prior,
                )
        elif precursor_realism_scorer is not None:
            batch = rank_precursor_realism(batch, policy=policy)
        for candidate in batch:
            if candidate.strategy_id:
                candidates.append(candidate)
            else:
                unresolved.setdefault(
                    (candidate.template_id, candidate.precursor_smiles),
                    candidate,
                )
        groups = group_strategy_candidates(
            candidates,
            top_k_strategies=max(1, len(candidates)),
            max_realizations_per_strategy=max_realizations_per_strategy,
        )
        strategies = groups[:top_k_strategies]
        if diversify:
            representatives = _apply_strategic_candidate_reserve(
                (group.representative for group in strategies),
                {
                    level: tuple(
                        group.representative
                        for group in groups
                        if group.representative.abstraction_level == level
                    )
                    for level, _ in level_diagnostics
                },
                top_k=top_k_strategies,
                policy=policy,
            )
            by_id = {group.strategy_id: group for group in groups}
            strategies = tuple(
                replace(
                    by_id[candidate.strategy_id],
                    strategy_rank=rank,
                    representative=candidate,
                )
                for rank, candidate in enumerate(representatives, start=1)
            )
        if len(strategies) >= top_k_strategies:
            break
    return StrategySearchResult(
        strategies=strategies,
        diagnostics=OperatorLadderDiagnostics(
            levels_attempted=tuple(level for level, _ in level_diagnostics),
            level_diagnostics=tuple(level_diagnostics),
        ),
        requested_strategy_count=top_k_strategies,
        max_templates_per_level=max_templates_to_apply,
        max_validations_per_level=max_candidates_to_validate,
        max_realizations_per_strategy=max_realizations_per_strategy,
        unresolved_candidates=tuple(unresolved.values()),
    )


def disconnect_strategies(
    target_smiles: str,
    library: GenericTemplateLibrary,
    *,
    top_k_strategies: int = 10,
    max_realizations_per_strategy: int = 3,
    max_templates_to_apply: int = 500,
    max_candidates_to_validate: int = 100,
    use_context: bool = True,
    include_l0: bool = True,
    diversify: bool = True,
    use_hierarchical_ranking: bool = True,
    precursor_realism_scorer: (
        Callable[[str], tuple[PrecursorRealismAssessment, ...]] | None
    ) = None,
) -> tuple[StrategyProposal, ...]:
    """Return verified strategies; use the detailed API for search diagnostics."""

    return disconnect_strategies_detailed(
        target_smiles,
        library,
        top_k_strategies=top_k_strategies,
        max_realizations_per_strategy=max_realizations_per_strategy,
        max_templates_to_apply=max_templates_to_apply,
        max_candidates_to_validate=max_candidates_to_validate,
        use_context=use_context,
        include_l0=include_l0,
        diversify=diversify,
        use_hierarchical_ranking=use_hierarchical_ranking,
        precursor_realism_scorer=precursor_realism_scorer,
    ).strategies


__all__ = [
    "StrategySearchResult",
    "disconnect_strategies",
    "disconnect_strategies_detailed",
    "group_strategy_candidates",
]
