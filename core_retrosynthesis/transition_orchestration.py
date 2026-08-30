"""Bounded orchestration over deterministic single-step transition providers.

The agent-facing action contains only stable identifiers and a proposal budget.
Molecular structures remain inside the deterministic environment, and every
provider result crosses the same admission boundary before it can be consumed
by a route search or shown as an actionable proposal.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Mapping, Protocol, Sequence

from reactive_taxonomy import featurize_reaction

from .chemistry import canonical_smiles, digest
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
)
from .generic_search import disconnect_operator_ladder_detailed


TRANSITION_ORCHESTRATION_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class TransitionProviderMetadata:
    """Stable identity and bounded capabilities of one transition provider."""

    provider_id: str
    display_name: str
    capability_tags: tuple[str, ...] = ()
    maximum_proposal_budget: int = 10

    def __post_init__(self) -> None:
        if not self.provider_id or not self.display_name:
            raise ValueError("transition provider identity must be nonempty")
        if self.maximum_proposal_budget < 1:
            raise ValueError("provider proposal budget must be positive")


@dataclass(frozen=True)
class TransitionProviderBatch:
    """Raw candidates and provider-local diagnostics from one bounded call."""

    candidates: tuple[GenericDisconnectionCandidate, ...]
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


class TransitionProvider(Protocol):
    """Narrow contract implemented by deterministic transition generators."""

    @property
    def metadata(self) -> TransitionProviderMetadata: ...

    def expand(
        self,
        target_smiles: str,
        proposal_budget: int,
    ) -> TransitionProviderBatch: ...


@dataclass(frozen=True)
class CallableTransitionProvider:
    """Adapt an existing deterministic expansion callable to the provider API."""

    metadata: TransitionProviderMetadata
    expansion_function: Callable[
        [str, int],
        tuple[GenericDisconnectionCandidate, ...] | TransitionProviderBatch,
    ]

    def expand(
        self,
        target_smiles: str,
        proposal_budget: int,
    ) -> TransitionProviderBatch:
        """Execute the wrapped function and normalize its bounded result."""

        result = self.expansion_function(target_smiles, proposal_budget)
        if isinstance(result, TransitionProviderBatch):
            return TransitionProviderBatch(
                candidates=result.candidates[:proposal_budget],
                diagnostics=result.diagnostics,
            )
        return TransitionProviderBatch(candidates=tuple(result[:proposal_budget]))


@dataclass(frozen=True)
class OperatorLadderTransitionProvider:
    """Canonical generic-operator ladder exposed as a transition provider."""

    library: GenericTemplateLibrary
    provider_metadata: TransitionProviderMetadata = field(
        default_factory=lambda: TransitionProviderMetadata(
            provider_id="generic_operator_ladder",
            display_name="Generic operator ladder",
            capability_tags=(
                "structure_backed",
                "forward_validated",
                "compatibility_assessed",
            ),
            maximum_proposal_budget=10,
        )
    )
    max_templates_to_apply: int = 300
    max_candidates_to_validate: int = 50
    use_context: bool = True
    include_l0: bool = True
    diversify: bool = True
    use_hierarchical_ranking: bool = True
    minimum_candidates_per_level: int = 0

    @property
    def metadata(self) -> TransitionProviderMetadata:
        """Return immutable provider capability metadata."""

        return self.provider_metadata

    def expand(
        self,
        target_smiles: str,
        proposal_budget: int,
    ) -> TransitionProviderBatch:
        """Run the existing chemistry-first ladder without bypassing validation."""

        candidates, diagnostics = disconnect_operator_ladder_detailed(
            target_smiles,
            self.library,
            top_k=proposal_budget,
            max_templates_to_apply=self.max_templates_to_apply,
            max_candidates_to_validate=self.max_candidates_to_validate,
            use_context=self.use_context,
            include_l0=self.include_l0,
            diversify=self.diversify,
            use_hierarchical_ranking=self.use_hierarchical_ranking,
            minimum_candidates_per_level=self.minimum_candidates_per_level,
            lazy_validation=True,
        )
        return TransitionProviderBatch(
            candidates=candidates,
            diagnostics=diagnostics.to_dict(),
        )


@dataclass(frozen=True)
class ExpansionLeaf:
    """Agent-visible leaf identity paired with environment-owned chemistry."""

    leaf_id: str
    molecule_smiles: str
    expandable: bool = True

    def __post_init__(self) -> None:
        canonical = canonical_smiles(self.molecule_smiles)
        if not self.leaf_id:
            raise ValueError("expansion leaf ID must be nonempty")
        if canonical is None or "." in canonical:
            raise ValueError("expansion leaf must contain one valid molecule")
        object.__setattr__(self, "molecule_smiles", canonical)


@dataclass(frozen=True)
class ExpansionState:
    """Small deterministic state projection accepted by the orchestrator."""

    state_id: str
    leaves: tuple[ExpansionLeaf, ...]

    def __post_init__(self) -> None:
        if not self.state_id:
            raise ValueError("expansion state ID must be nonempty")
        leaf_ids = tuple(leaf.leaf_id for leaf in self.leaves)
        if len(set(leaf_ids)) != len(leaf_ids):
            raise ValueError("expansion state leaf IDs must be unique")

    def leaf(self, leaf_id: str) -> ExpansionLeaf:
        """Resolve one stable leaf ID without exposing structure in the action."""

        for leaf in self.leaves:
            if leaf.leaf_id == leaf_id:
                return leaf
        raise ValueError(f"unknown expansion leaf ID: {leaf_id}")


@dataclass(frozen=True)
class ExpandLeafAction:
    """Bounded agent action containing identifiers rather than structures."""

    leaf_id: str
    provider_id: str
    proposal_budget: int
    action_id: str = ""

    def __post_init__(self) -> None:
        if not self.leaf_id or not self.provider_id:
            raise ValueError("expansion action IDs must be nonempty")
        if self.proposal_budget < 1:
            raise ValueError("proposal budget must be positive")
        expected = digest(
            "EXPANDACTION1",
            self.leaf_id,
            self.provider_id,
            str(self.proposal_budget),
        )
        if self.action_id and self.action_id != expected:
            raise ValueError("expansion action ID contradicts its fields")
        if not self.action_id:
            object.__setattr__(self, "action_id", expected)


@dataclass(frozen=True)
class RejectedTransition:
    """One provider result rejected at the deterministic admission boundary."""

    provider_id: str
    provider_rank: int
    reason: str
    precursor_smiles: str


@dataclass(frozen=True)
class AdmittedTransition:
    """One verified transition with explicit provider provenance."""

    transition_id: str
    provider_id: str
    provider_rank: int
    candidate: GenericDisconnectionCandidate

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible transition record."""

        return {
            "transition_id": self.transition_id,
            "provider_id": self.provider_id,
            "provider_rank": self.provider_rank,
            "candidate": self.candidate.to_dict(),
        }


@dataclass(frozen=True)
class ProviderExpansionOutcome:
    """Auditable result of one action after deterministic admission."""

    action: ExpandLeafAction
    target_smiles: str
    raw_candidate_count: int
    admitted: tuple[AdmittedTransition, ...]
    rejected: tuple[RejectedTransition, ...]
    duplicate_candidate_count: int
    provider_diagnostics: Mapping[str, Any] = field(default_factory=dict)
    schema_version: str = TRANSITION_ORCHESTRATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible provider outcome."""

        return {
            "action": asdict(self.action),
            "target_smiles": self.target_smiles,
            "raw_candidate_count": self.raw_candidate_count,
            "admitted": [item.to_dict() for item in self.admitted],
            "rejected": [asdict(item) for item in self.rejected],
            "duplicate_candidate_count": self.duplicate_candidate_count,
            "provider_diagnostics": dict(self.provider_diagnostics),
            "schema_version": self.schema_version,
        }


@dataclass(frozen=True)
class ShadowTransition:
    """Deduplicated transition and every provider that proposed it."""

    transition_id: str
    candidate: GenericDisconnectionCandidate
    provider_ranks: tuple[tuple[str, int], ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible cross-provider transition."""

        return {
            "transition_id": self.transition_id,
            "candidate": self.candidate.to_dict(),
            "provider_ranks": dict(self.provider_ranks),
        }


@dataclass(frozen=True)
class ShadowExpansionReport:
    """Read-only comparison of providers on the same route leaf."""

    state_id: str
    leaf_id: str
    outcomes: tuple[ProviderExpansionOutcome, ...]
    unique_transitions: tuple[ShadowTransition, ...]
    shared_transition_count: int
    schema_version: str = TRANSITION_ORCHESTRATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible shadow-mode report."""

        return {
            "state_id": self.state_id,
            "leaf_id": self.leaf_id,
            "outcomes": [outcome.to_dict() for outcome in self.outcomes],
            "unique_transitions": [
                transition.to_dict() for transition in self.unique_transitions
            ],
            "shared_transition_count": self.shared_transition_count,
            "schema_version": self.schema_version,
        }


def _transition_id(candidate: GenericDisconnectionCandidate) -> str:
    return digest(
        "TRANSITION1",
        candidate.target_smiles,
        candidate.precursor_smiles,
        candidate.operator_id,
        candidate.disconnection_site_key,
    )


def _candidate_rejection_reason(
    candidate: GenericDisconnectionCandidate,
    target_smiles: str,
) -> str | None:
    if candidate.forward_validation_status != "verified_signature":
        return "provider_candidate_not_signature_verified"
    candidate_target = canonical_smiles(candidate.target_smiles)
    if candidate_target != target_smiles:
        return "provider_candidate_target_mismatch"
    precursor = canonical_smiles(candidate.precursor_smiles)
    if precursor is None:
        return "provider_candidate_invalid_precursors"
    parts = candidate.proposed_reaction_smiles.split(">>")
    if len(parts) != 2:
        return "provider_candidate_invalid_reaction"
    reaction_precursors = canonical_smiles(parts[0])
    reaction_target = canonical_smiles(parts[1])
    if reaction_precursors != precursor or reaction_target != target_smiles:
        return "provider_candidate_reaction_mismatch"
    if candidate.precursor_compatibility_disposition == "reject":
        return "precursor_compatibility_rejected"
    if candidate.reaction_compatibility_disposition == "reject":
        return "reaction_compatibility_rejected"
    validation_reaction = (
        candidate.condition_query_reaction_smiles
        or candidate.proposed_reaction_smiles
    )
    validation_parts = validation_reaction.split(">>")
    if len(validation_parts) != 2:
        return "provider_candidate_invalid_validation_reaction"
    validation_precursors = canonical_smiles(validation_parts[0])
    validation_target = canonical_smiles(validation_parts[1])
    if validation_precursors != precursor or validation_target != target_smiles:
        return "provider_candidate_validation_reaction_mismatch"
    analysis = featurize_reaction(validation_reaction)
    if not analysis.valid or analysis.reaction_signature is None:
        return "provider_candidate_signature_not_reproduced"
    return None


class TransitionProviderOrchestrator:
    """Resolve bounded actions and compare providers without mutating a tree."""

    def __init__(self, providers: Sequence[TransitionProvider]) -> None:
        values = tuple(providers)
        provider_ids = tuple(provider.metadata.provider_id for provider in values)
        if not values:
            raise ValueError("at least one transition provider is required")
        if len(set(provider_ids)) != len(provider_ids):
            raise ValueError("transition provider IDs must be unique")
        self._providers = {
            provider.metadata.provider_id: provider for provider in values
        }

    @property
    def provider_metadata(self) -> tuple[TransitionProviderMetadata, ...]:
        """Return provider metadata in deterministic registration order."""

        return tuple(provider.metadata for provider in self._providers.values())

    def expand(
        self,
        state: ExpansionState,
        action: ExpandLeafAction,
    ) -> ProviderExpansionOutcome:
        """Execute one identifier-only action and admit valid transitions."""

        leaf = state.leaf(action.leaf_id)
        if not leaf.expandable:
            raise ValueError("selected expansion leaf is not expandable")
        provider = self._providers.get(action.provider_id)
        if provider is None:
            raise ValueError(f"unknown transition provider ID: {action.provider_id}")
        if action.proposal_budget > provider.metadata.maximum_proposal_budget:
            raise ValueError("action proposal budget exceeds provider capability")

        batch = provider.expand(leaf.molecule_smiles, action.proposal_budget)
        raw_candidates = tuple(batch.candidates[: action.proposal_budget])
        admitted: list[AdmittedTransition] = []
        rejected: list[RejectedTransition] = []
        seen: set[str] = set()
        duplicate_count = 0
        for rank, candidate in enumerate(raw_candidates, start=1):
            reason = _candidate_rejection_reason(candidate, leaf.molecule_smiles)
            if reason is not None:
                rejected.append(
                    RejectedTransition(
                        provider_id=action.provider_id,
                        provider_rank=rank,
                        reason=reason,
                        precursor_smiles=candidate.precursor_smiles,
                    )
                )
                continue
            transition_id = _transition_id(candidate)
            if transition_id in seen:
                duplicate_count += 1
                continue
            seen.add(transition_id)
            admitted.append(
                AdmittedTransition(
                    transition_id=transition_id,
                    provider_id=action.provider_id,
                    provider_rank=rank,
                    candidate=candidate,
                )
            )
        return ProviderExpansionOutcome(
            action=action,
            target_smiles=leaf.molecule_smiles,
            raw_candidate_count=len(raw_candidates),
            admitted=tuple(admitted),
            rejected=tuple(rejected),
            duplicate_candidate_count=duplicate_count,
            provider_diagnostics=batch.diagnostics,
        )

    def compare_shadow(
        self,
        state: ExpansionState,
        *,
        leaf_id: str,
        provider_ids: Sequence[str],
        proposal_budget: int,
    ) -> ShadowExpansionReport:
        """Call selected providers and return a non-mutating comparison report."""

        selected_ids = tuple(provider_ids)
        if not selected_ids:
            raise ValueError("shadow comparison requires at least one provider")
        if len(set(selected_ids)) != len(selected_ids):
            raise ValueError("shadow comparison provider IDs must be unique")
        outcomes = tuple(
            self.expand(
                state,
                ExpandLeafAction(
                    leaf_id=leaf_id,
                    provider_id=provider_id,
                    proposal_budget=proposal_budget,
                ),
            )
            for provider_id in selected_ids
        )
        merged: dict[str, tuple[GenericDisconnectionCandidate, list[tuple[str, int]]]] = {}
        for outcome in outcomes:
            for transition in outcome.admitted:
                current = merged.get(transition.transition_id)
                if current is None:
                    merged[transition.transition_id] = (
                        transition.candidate,
                        [(transition.provider_id, transition.provider_rank)],
                    )
                else:
                    current[1].append(
                        (transition.provider_id, transition.provider_rank)
                    )
        unique = tuple(
            ShadowTransition(
                transition_id=transition_id,
                candidate=candidate,
                provider_ranks=tuple(provider_ranks),
            )
            for transition_id, (candidate, provider_ranks) in merged.items()
        )
        return ShadowExpansionReport(
            state_id=state.state_id,
            leaf_id=leaf_id,
            outcomes=outcomes,
            unique_transitions=unique,
            shared_transition_count=sum(
                len(transition.provider_ranks) > 1 for transition in unique
            ),
        )


__all__ = [
    "TRANSITION_ORCHESTRATION_SCHEMA_VERSION",
    "AdmittedTransition",
    "CallableTransitionProvider",
    "ExpandLeafAction",
    "ExpansionLeaf",
    "ExpansionState",
    "OperatorLadderTransitionProvider",
    "ProviderExpansionOutcome",
    "RejectedTransition",
    "ShadowExpansionReport",
    "ShadowTransition",
    "TransitionProvider",
    "TransitionProviderBatch",
    "TransitionProviderMetadata",
    "TransitionProviderOrchestrator",
]
