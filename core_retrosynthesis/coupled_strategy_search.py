"""Exact two-step strategy replay for v1/v2 retrosynthesis comparison.

This module deliberately implements the smallest executable contract for the
coupled-route POC.  It replays an observed strategy only for its exact final
product, validates the producer-to-consumer intermediate, and always exposes
the two physical reactions.  It does not claim that an observed pair is yet a
transferable reaction template.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping, Sequence

from .chemistry import canonical_smiles


COUPLED_STRATEGY_SEARCH_SCHEMA_VERSION = "1.0"
COUPLED_STRATEGY_SEARCH_ALGORITHM_VERSION = "coupled_strategy_replay.v1"
CoupledStrategyPolicyVersion = Literal["v1", "v2"]
POLICY_VERSIONS = frozenset({"v1", "v2"})
V1_ADMITTED_RELATIONSHIPS = frozenset(
    {"handle_progression", "same_site_coupled"}
)
V2_ADMITTED_DEPENDENCIES = frozenset(
    {
        "created_handle_consumed",
        "activation_then_conversion",
        "continued_site_transformation",
    }
)


@dataclass(frozen=True)
class CoupledStrategyPhysicalStep:
    """One independently inspectable reaction inside a logical action."""

    forward_step_number: int
    reaction_smiles: str
    reactants_smiles: str
    product_smiles: str
    source_reaction_id: str | None

    def __post_init__(self) -> None:
        if self.forward_step_number not in {1, 2}:
            raise ValueError("physical step number must be one or two")
        if not self.reaction_smiles:
            raise ValueError("physical step requires a reaction")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible physical step."""

        return asdict(self)


@dataclass(frozen=True)
class CoupledStrategyReplayAction:
    """One exact-target logical action expanded into two physical reactions."""

    policy_version: CoupledStrategyPolicyVersion
    occurrence_id: str
    exact_strategy_id: str
    typed_strategy_id: str
    relationship_class: str
    dependency_class: str
    target_smiles: str
    intermediate_smiles: str
    terminal_precursor_smiles: str
    logical_reaction_smiles: str
    physical_steps: tuple[CoupledStrategyPhysicalStep, ...]
    coupling_score: float
    patent_id: str | None
    source_route_id: str | None
    validation_status: str = "validated_exact_replay"
    schema_version: str = COUPLED_STRATEGY_SEARCH_SCHEMA_VERSION
    algorithm_version: str = COUPLED_STRATEGY_SEARCH_ALGORITHM_VERSION

    def __post_init__(self) -> None:
        if self.policy_version not in POLICY_VERSIONS:
            raise ValueError("unsupported coupled-strategy policy")
        if len(self.physical_steps) != 2:
            raise ValueError("coupled strategy must retain two physical steps")
        if tuple(step.forward_step_number for step in self.physical_steps) != (
            1,
            2,
        ):
            raise ValueError("physical steps must remain in forward order")
        if not 0.0 <= self.coupling_score <= 1.0:
            raise ValueError("coupling score must be between zero and one")

    @property
    def retrosynthetic_steps(self) -> tuple[CoupledStrategyPhysicalStep, ...]:
        """Return the physical reactions in target-to-starting-material order."""

        return tuple(reversed(self.physical_steps))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible logical action."""

        return {
            "policy_version": self.policy_version,
            "occurrence_id": self.occurrence_id,
            "exact_strategy_id": self.exact_strategy_id,
            "typed_strategy_id": self.typed_strategy_id,
            "relationship_class": self.relationship_class,
            "dependency_class": self.dependency_class,
            "target_smiles": self.target_smiles,
            "intermediate_smiles": self.intermediate_smiles,
            "terminal_precursor_smiles": self.terminal_precursor_smiles,
            "logical_reaction_smiles": self.logical_reaction_smiles,
            "physical_steps": [step.to_dict() for step in self.physical_steps],
            "retrosynthetic_step_numbers": [
                step.forward_step_number for step in self.retrosynthetic_steps
            ],
            "coupling_score": self.coupling_score,
            "patent_id": self.patent_id,
            "source_route_id": self.source_route_id,
            "validation_status": self.validation_status,
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
        }


@dataclass(frozen=True)
class CoupledStrategySearchDiagnostics:
    """Transparent counters for one exact-replay search."""

    source_occurrence_count: int
    target_match_count: int
    policy_admitted_count: int
    policy_rejected_count: int
    invalid_chain_count: int
    valid_action_count: int
    returned_action_count: int

    def to_dict(self) -> dict[str, int]:
        """Return JSON-compatible search diagnostics."""

        return asdict(self)


@dataclass(frozen=True)
class CoupledStrategySearchResult:
    """Ranked exact-replay actions plus search diagnostics."""

    policy_version: CoupledStrategyPolicyVersion
    query_target_smiles: str
    actions: tuple[CoupledStrategyReplayAction, ...]
    diagnostics: CoupledStrategySearchDiagnostics
    schema_version: str = COUPLED_STRATEGY_SEARCH_SCHEMA_VERSION
    algorithm_version: str = COUPLED_STRATEGY_SEARCH_ALGORITHM_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible search result."""

        return {
            "policy_version": self.policy_version,
            "query_target_smiles": self.query_target_smiles,
            "actions": [action.to_dict() for action in self.actions],
            "diagnostics": self.diagnostics.to_dict(),
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
        }


@dataclass(frozen=True)
class CoupledStrategyPolicyComparison:
    """Paired v1/v2 exact-replay results for one target."""

    query_target_smiles: str
    v1: CoupledStrategySearchResult
    v2: CoupledStrategySearchResult

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible paired comparison."""

        return {
            "query_target_smiles": self.query_target_smiles,
            "v1": self.v1.to_dict(),
            "v2": self.v2.to_dict(),
        }


def _reaction_sides(reaction_smiles: str) -> tuple[str, str] | None:
    if reaction_smiles.count(">>") == 1:
        reactants, products = reaction_smiles.split(">>", 1)
    elif reaction_smiles.count(">") == 2:
        reactants, _, products = reaction_smiles.split(">", 2)
    else:
        return None
    if not reactants or not products:
        return None
    return reactants, products


def _canonical_components(smiles: str) -> tuple[str, ...] | None:
    values = []
    for component in smiles.split("."):
        canonical = canonical_smiles(component)
        if canonical is None:
            return None
        values.append(canonical)
    return tuple(values)


def _policy_admits(
    occurrence: Mapping[str, Any], policy_version: CoupledStrategyPolicyVersion
) -> bool:
    if policy_version == "v1":
        return (
            str(occurrence.get("relationship_class") or "")
            in V1_ADMITTED_RELATIONSHIPS
        )
    return (
        str(occurrence.get("dependency_class") or "")
        in V2_ADMITTED_DEPENDENCIES
        and str(occurrence.get("admission_class") or "") == "strict"
    )


def _validated_action(
    occurrence: Mapping[str, Any],
    *,
    query_target: str,
    policy_version: CoupledStrategyPolicyVersion,
) -> CoupledStrategyReplayAction | None:
    first_reaction = str(occurrence.get("first_reaction_smiles") or "")
    second_reaction = str(occurrence.get("second_reaction_smiles") or "")
    first_sides = _reaction_sides(first_reaction)
    second_sides = _reaction_sides(second_reaction)
    intermediate = canonical_smiles(
        str(occurrence.get("intermediate_smiles") or "")
    )
    target = canonical_smiles(str(occurrence.get("final_product_smiles") or ""))
    if (
        first_sides is None
        or second_sides is None
        or intermediate is None
        or target is None
        or target != query_target
    ):
        return None
    first_reactants, first_product = first_sides
    second_reactants, second_product = second_sides
    canonical_first_reactants = canonical_smiles(first_reactants)
    canonical_first_product = canonical_smiles(first_product)
    canonical_second_reactants = canonical_smiles(second_reactants)
    canonical_second_product = canonical_smiles(second_product)
    second_components = _canonical_components(second_reactants)
    if (
        canonical_first_reactants is None
        or canonical_second_reactants is None
        or canonical_first_product != intermediate
        or canonical_second_product != target
        or second_components is None
        or intermediate not in second_components
    ):
        return None
    remaining_second = list(second_components)
    remaining_second.remove(intermediate)
    terminal_components = [canonical_first_reactants, *remaining_second]
    terminal_precursors = canonical_smiles(".".join(terminal_components))
    if terminal_precursors is None:
        return None
    logical_reaction = f"{terminal_precursors}>>{target}"
    physical_steps = (
        CoupledStrategyPhysicalStep(
            forward_step_number=1,
            reaction_smiles=first_reaction,
            reactants_smiles=canonical_first_reactants,
            product_smiles=intermediate,
            source_reaction_id=(
                str(occurrence["first_source_reaction_id"])
                if occurrence.get("first_source_reaction_id") is not None
                else None
            ),
        ),
        CoupledStrategyPhysicalStep(
            forward_step_number=2,
            reaction_smiles=second_reaction,
            reactants_smiles=canonical_second_reactants,
            product_smiles=target,
            source_reaction_id=(
                str(occurrence["second_source_reaction_id"])
                if occurrence.get("second_source_reaction_id") is not None
                else None
            ),
        ),
    )
    return CoupledStrategyReplayAction(
        policy_version=policy_version,
        occurrence_id=str(occurrence.get("occurrence_id") or ""),
        exact_strategy_id=str(occurrence.get("exact_strategy_id") or ""),
        typed_strategy_id=str(occurrence.get("typed_strategy_id") or ""),
        relationship_class=str(
            occurrence.get("relationship_class") or "unresolved"
        ),
        dependency_class=str(
            occurrence.get("dependency_class")
            or occurrence.get("relationship_class")
            or "unresolved"
        ),
        target_smiles=target,
        intermediate_smiles=intermediate,
        terminal_precursor_smiles=terminal_precursors,
        logical_reaction_smiles=logical_reaction,
        physical_steps=physical_steps,
        coupling_score=float(occurrence.get("coupling_score") or 0.0),
        patent_id=(
            str(occurrence["patent_id"])
            if occurrence.get("patent_id") is not None
            else None
        ),
        source_route_id=(
            str(occurrence["source_route_id"])
            if occurrence.get("source_route_id") is not None
            else None
        ),
    )


def _occurrences(
    source: Mapping[str, Any] | Iterable[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    if isinstance(source, Mapping):
        sample = source.get("sample")
        if not isinstance(sample, Sequence):
            raise ValueError("coupled-strategy report requires a sample sequence")
        return tuple(item for item in sample if isinstance(item, Mapping))
    return tuple(source)


def search_coupled_strategy_replays(
    target_smiles: str,
    source: Mapping[str, Any] | Iterable[Mapping[str, Any]],
    *,
    policy_version: CoupledStrategyPolicyVersion,
    top_k: int = 10,
) -> CoupledStrategySearchResult:
    """Return exact-target, chain-validated two-step strategy replays."""

    if policy_version not in POLICY_VERSIONS:
        raise ValueError("policy version must be 'v1' or 'v2'")
    if top_k < 1:
        raise ValueError("top-k must be positive")
    query_target = canonical_smiles(target_smiles)
    if query_target is None:
        raise ValueError("target SMILES is invalid")
    occurrences = _occurrences(source)
    target_matches = 0
    admitted = 0
    rejected = 0
    invalid = 0
    actions = []
    for occurrence in occurrences:
        occurrence_target = canonical_smiles(
            str(occurrence.get("final_product_smiles") or "")
        )
        if occurrence_target != query_target:
            continue
        target_matches += 1
        if not _policy_admits(occurrence, policy_version):
            rejected += 1
            continue
        admitted += 1
        action = _validated_action(
            occurrence,
            query_target=query_target,
            policy_version=policy_version,
        )
        if action is None:
            invalid += 1
            continue
        actions.append(action)
    actions.sort(
        key=lambda item: (
            -item.coupling_score,
            item.typed_strategy_id,
            item.occurrence_id,
        )
    )
    returned = tuple(actions[:top_k])
    diagnostics = CoupledStrategySearchDiagnostics(
        source_occurrence_count=len(occurrences),
        target_match_count=target_matches,
        policy_admitted_count=admitted,
        policy_rejected_count=rejected,
        invalid_chain_count=invalid,
        valid_action_count=len(actions),
        returned_action_count=len(returned),
    )
    return CoupledStrategySearchResult(
        policy_version=policy_version,
        query_target_smiles=query_target,
        actions=returned,
        diagnostics=diagnostics,
    )


def compare_coupled_strategy_policies(
    target_smiles: str,
    source: Mapping[str, Any] | Iterable[Mapping[str, Any]],
    *,
    top_k: int = 10,
) -> CoupledStrategyPolicyComparison:
    """Run paired v1/v2 exact replay over one shared occurrence collection."""

    occurrences = _occurrences(source)
    v1 = search_coupled_strategy_replays(
        target_smiles,
        occurrences,
        policy_version="v1",
        top_k=top_k,
    )
    v2 = search_coupled_strategy_replays(
        target_smiles,
        occurrences,
        policy_version="v2",
        top_k=top_k,
    )
    return CoupledStrategyPolicyComparison(
        query_target_smiles=v1.query_target_smiles,
        v1=v1,
        v2=v2,
    )


def write_coupled_strategy_policy_comparison(
    comparison: CoupledStrategyPolicyComparison, output_path: str | Path
) -> dict[str, Any]:
    """Write a paired replay comparison and return a compact summary."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(comparison.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {
        "output_json": str(destination.resolve()),
        "json_bytes": destination.stat().st_size,
        "query_target_smiles": comparison.query_target_smiles,
        "v1_action_count": len(comparison.v1.actions),
        "v2_action_count": len(comparison.v2.actions),
    }


__all__ = [
    "COUPLED_STRATEGY_SEARCH_ALGORITHM_VERSION",
    "COUPLED_STRATEGY_SEARCH_SCHEMA_VERSION",
    "POLICY_VERSIONS",
    "V1_ADMITTED_RELATIONSHIPS",
    "V2_ADMITTED_DEPENDENCIES",
    "CoupledStrategyPhysicalStep",
    "CoupledStrategyPolicyComparison",
    "CoupledStrategyReplayAction",
    "CoupledStrategySearchDiagnostics",
    "CoupledStrategySearchResult",
    "compare_coupled_strategy_policies",
    "search_coupled_strategy_replays",
    "write_coupled_strategy_policy_comparison",
]
