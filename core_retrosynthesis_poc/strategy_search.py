"""Strategy-grouped single-step retrosynthesis over validated candidates."""

from __future__ import annotations

from typing import Callable, Iterable

from cas_tools import PrecursorRealismAssessment

from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    StrategyProposal,
)
from .generic_search import disconnect_operator_ladder
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
            raise ValueError(
                "strategy grouping requires verified-signature candidates"
            )
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
    """Return top-k distinct strategies with bounded concrete realizations.

    This first implementation groups only after the ordinary hard forward-graph
    and operator-signature validation path.  The larger flat pool prevents
    tactical variants from consuming the strategy result budget while leaving
    the established candidate generator and flat API unchanged.
    """

    if top_k_strategies < 1:
        raise ValueError("top-k strategies must be positive")
    if max_realizations_per_strategy < 1:
        raise ValueError("maximum realizations per strategy must be positive")
    if max_candidates_to_validate < 1:
        raise ValueError("maximum candidates to validate must be positive")

    policy = load_retrosynthesis_ranking_policy()
    flat_pool_size = min(
        max_candidates_to_validate,
        max(
            top_k_strategies,
            top_k_strategies * policy.candidate_pool_multiplier,
            top_k_strategies * max_realizations_per_strategy,
        ),
    )
    candidates = disconnect_operator_ladder(
        target_smiles,
        library,
        top_k=flat_pool_size,
        max_templates_to_apply=max_templates_to_apply,
        max_candidates_to_validate=max_candidates_to_validate,
        use_context=use_context,
        include_l0=include_l0,
        diversify=diversify,
        use_hierarchical_ranking=use_hierarchical_ranking,
        precursor_realism_scorer=precursor_realism_scorer,
    )
    return group_strategy_candidates(
        candidates,
        top_k_strategies=top_k_strategies,
        max_realizations_per_strategy=max_realizations_per_strategy,
    )


__all__ = ["disconnect_strategies", "group_strategy_candidates"]
