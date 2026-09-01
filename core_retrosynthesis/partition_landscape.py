"""Build role-neutral partition landscapes from validated one-step operators."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .chemistry import canonical_smiles
from .generic_models import GenericDisconnectionCandidate
from .partition_projection import project_reaction_to_target
from .synthetic_partition import (
    InterfaceHypothesis,
    PartitionSearchDiagnostics,
    SyntheticPartition,
    SyntheticPartitionLandscape,
    analyze_partition_target,
    create_synthetic_partition,
)


SYNTHETIC_PARTITION_POLICY_PATH = (
    Path(__file__).with_name("definitions") / "synthetic_partition_policy.v1.json"
)


@dataclass(frozen=True)
class SyntheticPartitionPolicy:
    """Validated deterministic bounds and ordering for landscape generation."""

    definition_id: str
    schema_version: str
    minimum_k: int
    maximum_k: int
    maximum_seed_candidates: int
    maximum_combination_size: int
    maximum_combinations: int
    maximum_partitions_per_k: int
    maximum_returned_partitions: int
    close_context_threshold: float
    strong_independent_reference_support: int
    evidence_order: tuple[str, ...]
    source_order: tuple[str, ...]


@dataclass(frozen=True)
class _CandidateSeed:
    module_atom_sets: tuple[tuple[int, ...], ...]
    target_bonds: tuple[tuple[int, int, str], ...]
    operator_ids: tuple[str, ...]
    strategy_ids: tuple[str, ...]
    site_keys: tuple[str, ...]
    precedent_ids: tuple[str, ...]
    evidence_tokens: tuple[str, ...]
    evidence_level: str
    heuristic_score: float
    warnings: tuple[str, ...]

    @property
    def key(self) -> tuple[Any, ...]:
        return self.module_atom_sets, self.target_bonds

    def to_hypothesis(self) -> InterfaceHypothesis:
        """Return interface evidence without claiming joint realization."""

        return InterfaceHypothesis(
            interface_kind=(
                "validated_operator_boundary"
                if self.target_bonds
                else "validated_unary_state"
            ),
            target_bonds=self.target_bonds,
            candidate_operator_ids=self.operator_ids,
            candidate_strategy_ids=self.strategy_ids,
            precedent_reaction_ids=self.precedent_ids,
            evidence_level=self.evidence_level,
            warnings=self.warnings,
        )


@dataclass(frozen=True)
class _PartitionProposal:
    module_atom_sets: tuple[tuple[int, ...], ...]
    seeds: tuple[_CandidateSeed, ...]
    source_kind: str
    evidence_level: str
    heuristic_score: float
    warnings: tuple[str, ...]


def validate_synthetic_partition_policy(value: Mapping[str, Any]) -> list[str]:
    """Return stable validation errors for a partition policy definition."""

    errors: list[str] = []
    if value.get("definition_id") != "synthetic_partition_policy.v1":
        errors.append("invalid_definition_id")
    if value.get("schema_version") != "1.0":
        errors.append("invalid_schema_version")
    k_range = value.get("default_k_range")
    limits = value.get("search_limits")
    evidence = value.get("evidence")
    ranking = value.get("ranking")
    if not isinstance(k_range, Mapping):
        errors.append("missing_default_k_range")
        k_range = {}
    if not isinstance(limits, Mapping):
        errors.append("missing_search_limits")
        limits = {}
    if not isinstance(evidence, Mapping):
        errors.append("missing_evidence")
        evidence = {}
    if not isinstance(ranking, Mapping):
        errors.append("missing_ranking")
        ranking = {}
    try:
        minimum_k = int(k_range.get("minimum", 0))
        maximum_k = int(k_range.get("maximum", 0))
        if minimum_k < 1 or maximum_k < minimum_k:
            errors.append("invalid_k_range")
    except (TypeError, ValueError):
        errors.append("invalid_k_range")
    for field in (
        "maximum_seed_candidates",
        "maximum_combination_size",
        "maximum_combinations",
        "maximum_partitions_per_k",
        "maximum_returned_partitions",
    ):
        try:
            if int(limits.get(field, 0)) < 1:
                errors.append(f"invalid_{field}")
        except (TypeError, ValueError):
            errors.append(f"invalid_{field}")
    try:
        threshold = float(evidence.get("close_context_threshold", -1.0))
        if not 0.0 <= threshold <= 1.0:
            errors.append("invalid_close_context_threshold")
    except (TypeError, ValueError):
        errors.append("invalid_close_context_threshold")
    try:
        if int(evidence.get("strong_independent_reference_support", 0)) < 1:
            errors.append("invalid_strong_reference_support")
    except (TypeError, ValueError):
        errors.append("invalid_strong_reference_support")
    if tuple(ranking.get("evidence_order") or ()) != (
        "E4",
        "E3",
        "E2",
        "E1",
        "E0",
    ):
        errors.append("invalid_evidence_order")
    if not ranking.get("source_order"):
        errors.append("missing_source_order")
    return errors


@lru_cache(maxsize=1)
def load_synthetic_partition_policy() -> SyntheticPartitionPolicy:
    """Load and validate the versioned partition policy."""

    value = json.loads(SYNTHETIC_PARTITION_POLICY_PATH.read_text(encoding="utf-8"))
    errors = validate_synthetic_partition_policy(value)
    if errors:
        raise ValueError("invalid synthetic partition policy: " + ", ".join(errors))
    k_range = value["default_k_range"]
    limits = value["search_limits"]
    evidence = value["evidence"]
    ranking = value["ranking"]
    return SyntheticPartitionPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        minimum_k=int(k_range["minimum"]),
        maximum_k=int(k_range["maximum"]),
        maximum_seed_candidates=int(limits["maximum_seed_candidates"]),
        maximum_combination_size=int(limits["maximum_combination_size"]),
        maximum_combinations=int(limits["maximum_combinations"]),
        maximum_partitions_per_k=int(limits["maximum_partitions_per_k"]),
        maximum_returned_partitions=int(limits["maximum_returned_partitions"]),
        close_context_threshold=float(evidence["close_context_threshold"]),
        strong_independent_reference_support=int(
            evidence["strong_independent_reference_support"]
        ),
        evidence_order=tuple(str(item) for item in ranking["evidence_order"]),
        source_order=tuple(str(item) for item in ranking["source_order"]),
    )


def _candidate_evidence_level(
    candidate: GenericDisconnectionCandidate,
    policy: SyntheticPartitionPolicy,
) -> str:
    if (
        candidate.independent_reference_support
        >= policy.strong_independent_reference_support
        and candidate.context_similarity >= policy.close_context_threshold
    ):
        return "E3"
    return "E2"


def _stronger_evidence(
    left: str,
    right: str,
    policy: SyntheticPartitionPolicy,
) -> str:
    rank = {value: index for index, value in enumerate(policy.evidence_order)}
    return left if rank[left] <= rank[right] else right


def _merge_seed(
    left: _CandidateSeed,
    right: _CandidateSeed,
    policy: SyntheticPartitionPolicy,
) -> _CandidateSeed:
    return _CandidateSeed(
        module_atom_sets=left.module_atom_sets,
        target_bonds=left.target_bonds,
        operator_ids=tuple(sorted(set(left.operator_ids + right.operator_ids))),
        strategy_ids=tuple(sorted(set(left.strategy_ids + right.strategy_ids))),
        site_keys=tuple(sorted(set(left.site_keys + right.site_keys))),
        precedent_ids=tuple(sorted(set(left.precedent_ids + right.precedent_ids))),
        evidence_tokens=tuple(
            sorted(set(left.evidence_tokens + right.evidence_tokens))
        ),
        evidence_level=_stronger_evidence(
            left.evidence_level,
            right.evidence_level,
            policy,
        ),
        heuristic_score=max(left.heuristic_score, right.heuristic_score),
        warnings=tuple(sorted(set(left.warnings + right.warnings))),
    )


def _candidate_seed(
    target_smiles: str,
    candidate: GenericDisconnectionCandidate,
    policy: SyntheticPartitionPolicy,
) -> _CandidateSeed:
    if candidate.forward_validation_status != "verified_signature":
        raise ValueError("candidate_not_forward_verified")
    target = canonical_smiles(candidate.target_smiles)
    if target != target_smiles:
        raise ValueError("candidate_target_mismatch")
    reaction = (
        candidate.condition_query_reaction_smiles
        or candidate.proposed_reaction_smiles
    )
    projection = project_reaction_to_target(target_smiles, reaction)
    evidence_token = (
        candidate.strategy_id
        or candidate.realization_id
        or candidate.template_id
    )
    warnings = list(projection.warnings)
    if projection.mapping_evidence != "supplied_atom_mapping":
        warnings.append("ATOM_MAPPING_INFERRED_FOR_PARTITION_PROJECTION")
    return _CandidateSeed(
        module_atom_sets=projection.module_atom_sets,
        target_bonds=projection.target_boundary_bonds,
        operator_ids=tuple(item for item in (candidate.operator_id,) if item),
        strategy_ids=tuple(item for item in (candidate.strategy_id,) if item),
        site_keys=tuple(
            item for item in (candidate.disconnection_site_key,) if item
        ),
        precedent_ids=tuple(sorted(set(candidate.precedent_reaction_ids))),
        evidence_tokens=(evidence_token,),
        evidence_level=_candidate_evidence_level(candidate, policy),
        heuristic_score=round(float(candidate.score), 8),
        warnings=tuple(sorted(set(warnings))),
    )


def _refine_module_sets(
    target_atom_maps: tuple[int, ...],
    seeds: Sequence[_CandidateSeed],
) -> tuple[tuple[int, ...], ...]:
    signatures: dict[tuple[int, ...], list[int]] = {}
    for atom_map in target_atom_maps:
        signature = []
        for seed in seeds:
            memberships = tuple(
                index
                for index, block in enumerate(seed.module_atom_sets)
                if atom_map in block
            )
            if len(memberships) != 1:
                raise ValueError("candidate partition has invalid atom ownership")
            signature.append(memberships[0])
        signatures.setdefault(tuple(signature), []).append(atom_map)
    return tuple(
        sorted(
            (tuple(sorted(block)) for block in signatures.values()),
            key=lambda item: (item[0], len(item), item),
        )
    )


def _compatible_seed_combination(seeds: Sequence[_CandidateSeed]) -> bool:
    seen_bonds: set[tuple[int, int, str]] = set()
    seen_sites: set[str] = set()
    for seed in seeds:
        bonds = set(seed.target_bonds)
        sites = set(seed.site_keys)
        if bonds and seen_bonds.intersection(bonds):
            return False
        if sites and seen_sites.intersection(sites):
            return False
        seen_bonds.update(bonds)
        seen_sites.update(sites)
    return True


def _proposal_sort_key(
    proposal: _PartitionProposal,
    policy: SyntheticPartitionPolicy,
) -> tuple[Any, ...]:
    evidence_rank = {
        value: index for index, value in enumerate(policy.evidence_order)
    }
    source_rank = {value: index for index, value in enumerate(policy.source_order)}
    return (
        evidence_rank.get(proposal.evidence_level, len(evidence_rank)),
        source_rank.get(proposal.source_kind, len(source_rank)),
        -proposal.heuristic_score,
        len(proposal.module_atom_sets),
        proposal.module_atom_sets,
    )


def _partition_sort_key(
    partition: SyntheticPartition,
    policy: SyntheticPartitionPolicy,
) -> tuple[Any, ...]:
    evidence_rank = {
        value: index for index, value in enumerate(policy.evidence_order)
    }
    source_rank = {value: index for index, value in enumerate(policy.source_order)}
    return (
        evidence_rank.get(partition.evidence_level, len(evidence_rank)),
        source_rank.get(partition.source_kind, len(source_rank)),
        -partition.heuristic_score,
        partition.k,
        partition.partition_id,
    )


def build_operator_partition_landscape(
    target_smiles: str,
    candidates: Iterable[GenericDisconnectionCandidate],
    *,
    minimum_k: int | None = None,
    maximum_k: int | None = None,
    policy: SyntheticPartitionPolicy | None = None,
) -> SyntheticPartitionLandscape:
    """Build a bounded review-only landscape from validated operator results."""

    selected_policy = policy or load_synthetic_partition_policy()
    canonical, target_id, target_atoms, _ = analyze_partition_target(target_smiles)
    lower = selected_policy.minimum_k if minimum_k is None else int(minimum_k)
    upper = selected_policy.maximum_k if maximum_k is None else int(maximum_k)
    if lower < 1 or upper < lower:
        raise ValueError("invalid partition k range")
    supplied = tuple(candidates)
    rejection_counts: Counter[str] = Counter()
    rejected_candidate_count = 0
    seeds_by_key: dict[tuple[Any, ...], _CandidateSeed] = {}
    ordered_candidates = tuple(
        sorted(
            supplied,
            key=lambda item: (-item.score, item.template_id, item.precursor_smiles),
        )
    )[: selected_policy.maximum_seed_candidates]
    for candidate in ordered_candidates:
        try:
            seed = _candidate_seed(canonical, candidate, selected_policy)
        except ValueError as exc:
            rejection_counts[str(exc)] += 1
            rejected_candidate_count += 1
            continue
        if len(seed.module_atom_sets) > upper:
            rejection_counts["candidate_exceeds_maximum_k"] += 1
            rejected_candidate_count += 1
            continue
        current = seeds_by_key.get(seed.key)
        seeds_by_key[seed.key] = (
            seed
            if current is None
            else _merge_seed(current, seed, selected_policy)
        )
    seeds = tuple(
        sorted(
            seeds_by_key.values(),
            key=lambda item: (
                -item.heuristic_score,
                item.module_atom_sets,
                item.target_bonds,
            ),
        )
    )
    target_maps = tuple(atom.atom_map for atom in target_atoms)
    proposals: list[_PartitionProposal] = [
        _PartitionProposal(
            module_atom_sets=(target_maps,),
            seeds=(),
            source_kind="structural_baseline",
            evidence_level="E0",
            heuristic_score=0.0,
            warnings=("NO_DISCONNECTION_REQUIRED_OR_ESTABLISHED",),
        )
    ]
    for seed in seeds:
        if lower <= len(seed.module_atom_sets) <= upper or len(seed.module_atom_sets) == 1:
            proposals.append(
                _PartitionProposal(
                    module_atom_sets=seed.module_atom_sets,
                    seeds=(seed,),
                    source_kind="validated_operator_projection",
                    evidence_level=seed.evidence_level,
                    heuristic_score=seed.heuristic_score,
                    warnings=seed.warnings,
                )
            )
    generated_combinations = 0
    maximum_size = min(selected_policy.maximum_combination_size, len(seeds))
    for size in range(2, maximum_size + 1):
        for selected in combinations(seeds, size):
            if generated_combinations >= selected_policy.maximum_combinations:
                break
            generated_combinations += 1
            if not _compatible_seed_combination(selected):
                rejection_counts["incompatible_interface_combination"] += 1
                continue
            refined = _refine_module_sets(target_maps, selected)
            if len(refined) <= max(len(seed.module_atom_sets) for seed in selected):
                rejection_counts["combination_adds_no_partition_information"] += 1
                continue
            if not lower <= len(refined) <= upper:
                rejection_counts["combination_outside_k_range"] += 1
                continue
            proposals.append(
                _PartitionProposal(
                    module_atom_sets=refined,
                    seeds=tuple(selected),
                    source_kind="operator_combination_unrealized",
                    evidence_level="E1",
                    heuristic_score=round(
                        min(seed.heuristic_score for seed in selected),
                        8,
                    ),
                    warnings=("INTERFACE_COMBINATION_NOT_JOINTLY_REALIZED",),
                )
            )
        if generated_combinations >= selected_policy.maximum_combinations:
            break
    proposals_by_blocks: dict[
        tuple[tuple[int, ...], ...], list[_PartitionProposal]
    ] = {}
    for proposal in proposals:
        proposals_by_blocks.setdefault(proposal.module_atom_sets, []).append(proposal)
    partitions = []
    for blocks, grouped in proposals_by_blocks.items():
        ordered_group = tuple(
            sorted(grouped, key=lambda item: _proposal_sort_key(item, selected_policy))
        )
        best = ordered_group[0]
        all_seeds = {
            seed.key: seed for proposal in ordered_group for seed in proposal.seeds
        }
        warnings = set(best.warnings)
        if any(
            proposal.source_kind == "operator_combination_unrealized"
            for proposal in ordered_group[1:]
        ):
            warnings.add("ALTERNATIVE_INTERFACE_COMBINATION_UNREALIZED")
        partitions.append(
            create_synthetic_partition(
                canonical,
                blocks,
                source_kind=best.source_kind,
                evidence_level=best.evidence_level,
                interface_hypotheses=tuple(
                    seed.to_hypothesis()
                    for seed in sorted(
                        all_seeds.values(),
                        key=lambda item: (item.module_atom_sets, item.target_bonds),
                    )
                ),
                generation_evidence=(
                    token
                    for seed in all_seeds.values()
                    for token in seed.evidence_tokens
                ),
                heuristic_score=best.heuristic_score,
                warnings=warnings,
            )
        )
    ordered_partitions = tuple(
        sorted(partitions, key=lambda item: _partition_sort_key(item, selected_policy))
    )
    per_k: Counter[int] = Counter()
    selected_partitions = []
    for partition in ordered_partitions:
        if per_k[partition.k] >= selected_policy.maximum_partitions_per_k:
            continue
        if len(selected_partitions) >= selected_policy.maximum_returned_partitions:
            break
        selected_partitions.append(partition)
        per_k[partition.k] += 1
    partition_changing = tuple(seed for seed in seeds if len(seed.module_atom_sets) > 1)
    abstained = not partition_changing
    abstention_reasons = (
        ("NO_VALIDATED_PARTITION_CHANGING_INTERFACE",) if abstained else ()
    )
    coverage = ["structural_k1_baseline"]
    if seeds:
        coverage.append("validated_operator_projection")
    if any(len(proposal.seeds) > 1 for proposal in proposals):
        coverage.append("unrealized_interface_combinations")
    diagnostics = PartitionSearchDiagnostics(
        supplied_candidate_count=len(supplied),
        accepted_seed_count=len(seeds),
        rejected_candidate_count=rejected_candidate_count,
        generated_combination_count=generated_combinations,
        duplicate_partition_count=len(proposals) - len(proposals_by_blocks),
        returned_partition_count=len(selected_partitions),
        rejection_counts=tuple(sorted(rejection_counts.items())),
    )
    warnings = (
        "Operator combinations are partition hypotheses, not jointly realized routes.",
        "Heuristic scores are ordering aids, not calibrated probabilities.",
    )
    return SyntheticPartitionLandscape(
        target_id=target_id,
        target_smiles=canonical,
        target_atoms=target_atoms,
        partitions=tuple(selected_partitions),
        searched_k_values=tuple(sorted(set(range(lower, upper + 1)) | {1})),
        generation_coverage=tuple(coverage),
        unresolved_motifs=(),
        abstained=abstained,
        abstention_reasons=abstention_reasons,
        diagnostics=diagnostics,
        warnings=warnings,
    )


__all__ = [
    "SYNTHETIC_PARTITION_POLICY_PATH",
    "SyntheticPartitionPolicy",
    "build_operator_partition_landscape",
    "load_synthetic_partition_policy",
    "validate_synthetic_partition_policy",
]
