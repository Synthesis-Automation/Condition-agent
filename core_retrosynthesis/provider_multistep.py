"""Adapters between transition providers and canonical multistep route search."""

from __future__ import annotations

from typing import Optional

from .chemistry import canonical_smiles, digest
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
    OperatorLadderDiagnostics,
)
from .multistep import (
    MultistepRetrosynthesisRoute,
    OneStepExpansionBatch,
)
from .transition_orchestration import (
    AdmittedTransition,
    ExpandLeafAction,
    ExpansionLeaf,
    ExpansionState,
    ProviderExpansionOutcome,
    TransitionProviderOrchestrator,
)


def _operator_ladder_diagnostics(
    outcome: ProviderExpansionOutcome,
) -> Optional[OperatorLadderDiagnostics]:
    raw = outcome.provider_diagnostics
    raw_levels = raw.get("levels_attempted")
    raw_diagnostics = raw.get("level_diagnostics")
    if not isinstance(raw_levels, (list, tuple)) or not isinstance(
        raw_diagnostics, dict
    ):
        return None
    values = []
    try:
        for level in raw_levels:
            level_name = str(level)
            level_value = raw_diagnostics[level_name]
            if not isinstance(level_value, dict):
                return None
            values.append(
                (level_name, GenericSearchDiagnostics(**level_value))
            )
    except (KeyError, TypeError, ValueError):
        return None
    return OperatorLadderDiagnostics(
        levels_attempted=tuple(str(level) for level in raw_levels),
        level_diagnostics=tuple(values),
    )


class ProviderBackedOneStepExpander:
    """Expose one registered provider through the planner's expander contract.

    The adapter owns only call-local telemetry. It does not rank candidates or
    alter route chemistry; every returned candidate has already crossed the
    transition-provider admission boundary.
    """

    def __init__(
        self,
        orchestrator: TransitionProviderOrchestrator,
        provider_id: str,
    ) -> None:
        provider_ids = {
            item.provider_id for item in orchestrator.provider_metadata
        }
        if provider_id not in provider_ids:
            raise ValueError(f"unknown transition provider ID: {provider_id}")
        self._orchestrator = orchestrator
        self.provider_id = provider_id
        self._outcomes: dict[str, ProviderExpansionOutcome] = {}

    def __call__(
        self,
        product_smiles: str,
        top_k: int,
    ) -> OneStepExpansionBatch:
        """Expand one canonical product with a bounded identifier-only action."""

        canonical = canonical_smiles(product_smiles)
        if canonical is None or "." in canonical:
            raise ValueError("provider-backed expansion requires one valid product")
        leaf_id = digest("PROVIDERLEAF1", canonical)
        state = ExpansionState(
            state_id=digest("PROVIDERSTATE1", canonical),
            leaves=(ExpansionLeaf(leaf_id, canonical),),
        )
        outcome = self._orchestrator.expand(
            state,
            ExpandLeafAction(
                leaf_id=leaf_id,
                provider_id=self.provider_id,
                proposal_budget=top_k,
            ),
        )
        self._outcomes[canonical] = outcome
        return OneStepExpansionBatch(
            candidates=tuple(item.candidate for item in outcome.admitted),
            diagnostics=_operator_ladder_diagnostics(outcome),
        )

    @property
    def outcomes(self) -> tuple[ProviderExpansionOutcome, ...]:
        """Return outcomes in canonical target order for deterministic reports."""

        return tuple(self._outcomes[key] for key in sorted(self._outcomes))

    def outcome_for_target(
        self,
        target_smiles: str,
    ) -> Optional[ProviderExpansionOutcome]:
        """Return the cached provider outcome for one canonical target."""

        canonical = canonical_smiles(target_smiles)
        return self._outcomes.get(canonical or "")

    def attribution_for_candidate(
        self,
        candidate: GenericDisconnectionCandidate,
    ) -> Optional[AdmittedTransition]:
        """Resolve provider provenance for one candidate retained by a route."""

        outcome = self.outcome_for_target(candidate.target_smiles)
        if outcome is None:
            return None
        return next(
            (
                item
                for item in outcome.admitted
                if item.candidate.proposed_reaction_smiles
                == candidate.proposed_reaction_smiles
                and item.candidate.operator_id == candidate.operator_id
            ),
            None,
        )


def expansion_state_from_route(
    route: MultistepRetrosynthesisRoute,
) -> ExpansionState:
    """Project route leaves to stable-ID agent observations without new chemistry."""

    return ExpansionState(
        state_id=digest("ROUTEEXPANSIONSTATE1", route.route_id),
        leaves=tuple(
            ExpansionLeaf(
                leaf_id=leaf.route_node_id,
                molecule_smiles=leaf.canonical_smiles,
                expandable=(
                    not leaf.terminal and leaf.unresolved_reason is None
                ),
            )
            for leaf in route.leaves
        ),
    )


__all__ = [
    "ProviderBackedOneStepExpander",
    "expansion_state_from_route",
]

