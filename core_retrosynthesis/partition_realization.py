"""Bounded realization of role-neutral synthetic partitions.

This module guides the canonical multistep planner toward a selected target
partition.  It never creates reaction actions: every route edge still comes
from the existing validated one-step operator ladder.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from rdkit import Chem

from .chemistry import digest
from .generic_models import GenericTemplateLibrary
from .multistep import (
    ConditionEvidenceEvaluator,
    LiteratureLookup,
    MultistepGuidanceState,
    MultistepRetrosynthesisRoute,
    OneStepExpander,
    PrecursorRealismScorer,
    plan_multistep_routes,
)
from .partition_projection import (
    RouteFrontierPartition,
    RoutePartitionProjection,
    project_route_partitions,
)
from .route_action_policy import RouteActionPolicyModel
from .route_contract import ReactionRouteTree
from .route_refinement import RouteCandidateExclusion
from .route_tree import build_canonical_route_tree
from .route_verification import RouteVerificationReport, verify_planned_route
from .synthetic_partition import (
    PARTITION_EVIDENCE_LEVELS,
    PARTITION_REALIZATION_STATUSES,
    LatentModuleState,
    StrategicInterface,
    SyntheticPartition,
    analyze_partition_target,
    validate_synthetic_partition,
)


PARTITION_REALIZATION_SCHEMA_VERSION = "1.0"
PARTITION_REALIZATION_ALGORITHM_VERSION = "partition_guided_multistep.v1"
PARTITION_REALIZATION_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "synthetic_partition_realization_policy.v1.json"
)
_STATUS_ORDER = (
    "fully_realized",
    "partially_realized",
    "unrealized_but_plausible",
    "contradicted",
)
_EVIDENCE_ORDER = ("E4", "E3", "E2", "E1", "E0")


@dataclass(frozen=True)
class PartitionRealizationPolicy:
    """Validated defaults for bounded partition realization."""

    definition_id: str
    schema_version: str
    max_depth: int
    molecular_weight_threshold: float
    maximum_realizations: int
    route_candidate_limit: int
    per_step_top_k: int
    beam_width: int
    maximum_expansions: int
    maximum_templates_to_apply: int
    maximum_candidates_to_validate: int
    partial_interface_coverage_threshold: float
    status_order: tuple[str, ...]


@dataclass(frozen=True)
class PartitionRealizationConfig:
    """Concrete deterministic budget for one realization request."""

    max_depth: int = 3
    molecular_weight_threshold: float = 150.0
    maximum_realizations: int = 5
    route_candidate_limit: int = 20
    per_step_top_k: int = 8
    beam_width: int = 40
    maximum_expansions: int = 200
    maximum_templates_to_apply: int = 300
    maximum_candidates_to_validate: int = 50
    use_context: bool = True
    include_l0: bool = True
    diversify: bool = True
    use_hierarchical_ranking: bool = True

    def __post_init__(self) -> None:
        if not 1 <= self.max_depth <= 6:
            raise ValueError("partition realization depth must be in [1, 6]")
        if self.molecular_weight_threshold <= 0.0:
            raise ValueError("molecular-weight threshold must be positive")
        for value, name in (
            (self.maximum_realizations, "maximum realizations"),
            (self.route_candidate_limit, "route candidate limit"),
            (self.per_step_top_k, "per-step top-k"),
            (self.beam_width, "beam width"),
            (self.maximum_expansions, "maximum expansions"),
            (self.maximum_templates_to_apply, "maximum templates"),
            (
                self.maximum_candidates_to_validate,
                "maximum candidates to validate",
            ),
        ):
            if value < 1:
                raise ValueError(f"{name} must be positive")
        if self.route_candidate_limit < self.maximum_realizations:
            raise ValueError("route candidate limit must cover maximum realizations")
        if self.maximum_candidates_to_validate < self.per_step_top_k:
            raise ValueError("candidate validation budget must cover per-step top-k")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible search configuration."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionRealizationConfig":
        """Reconstruct a search configuration."""

        return cls(**{field: value[field] for field in cls.__dataclass_fields__})


@dataclass(frozen=True)
class PartitionFrontierMatch:
    """Permutation-invariant agreement between a route frontier and a partition."""

    frontier_id: str
    depth: int
    desired_k: int
    frontier_k: int
    exact_partition_match: bool
    matched_module_count: int
    atom_weighted_module_overlap: float
    boundary_precision: float
    boundary_recall: float
    boundary_f1: float
    desired_boundary_count: int
    matched_boundary_count: int
    extra_boundary_count: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible frontier match."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionFrontierMatch":
        """Reconstruct a frontier match."""

        return cls(**{field: value[field] for field in cls.__dataclass_fields__})


@dataclass(frozen=True)
class InterfaceRealization:
    """Route evidence for one requested strategic interface."""

    interface_id: str
    status: str
    realization_kind: str
    target_bonds: tuple[tuple[int, int, str], ...]
    matched_target_bonds: tuple[tuple[int, int, str], ...]
    route_frontier_ids: tuple[str, ...]
    route_step_ids: tuple[str, ...]
    operator_ids: tuple[str, ...]
    strategy_ids: tuple[str, ...]
    evidence_level: str
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.status not in {"realized", "unresolved", "contradicted"}:
            raise ValueError("unsupported interface realization status")
        if self.evidence_level not in PARTITION_EVIDENCE_LEVELS:
            raise ValueError("unsupported interface realization evidence")

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible interface evidence."""

        value = asdict(self)
        value["target_bonds"] = [list(item) for item in self.target_bonds]
        value["matched_target_bonds"] = [
            list(item) for item in self.matched_target_bonds
        ]
        value["route_frontier_ids"] = list(self.route_frontier_ids)
        value["route_step_ids"] = list(self.route_step_ids)
        value["operator_ids"] = list(self.operator_ids)
        value["strategy_ids"] = list(self.strategy_ids)
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "InterfaceRealization":
        """Reconstruct interface realization evidence."""

        return cls(
            interface_id=str(value["interface_id"]),
            status=str(value["status"]),
            realization_kind=str(value["realization_kind"]),
            target_bonds=tuple(
                (int(item[0]), int(item[1]), str(item[2]))
                for item in value.get("target_bonds") or ()
            ),
            matched_target_bonds=tuple(
                (int(item[0]), int(item[1]), str(item[2]))
                for item in value.get("matched_target_bonds") or ()
            ),
            route_frontier_ids=tuple(
                str(item) for item in value.get("route_frontier_ids") or ()
            ),
            route_step_ids=tuple(
                str(item) for item in value.get("route_step_ids") or ()
            ),
            operator_ids=tuple(str(item) for item in value.get("operator_ids") or ()),
            strategy_ids=tuple(str(item) for item in value.get("strategy_ids") or ()),
            evidence_level=str(value["evidence_level"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class RealizationEvidenceSummary:
    """Auditable structural evidence for one partition-guided route."""

    requested_interface_count: int
    realized_interface_count: int
    validated_interface_coverage: float
    exact_partition_match: bool
    atom_weighted_module_overlap: float
    boundary_precision: float
    boundary_recall: float
    boundary_f1: float
    weakest_interface_evidence: str | None
    validated_step_count: int
    unsupported_step_count: int
    dependency_graph_valid: bool
    forward_order_exists: bool
    terminal_leaf_count: int
    unresolved_leaf_count: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible evidence summary."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RealizationEvidenceSummary":
        """Reconstruct a realization evidence summary."""

        return cls(**{field: value[field] for field in cls.__dataclass_fields__})


@dataclass(frozen=True)
class StrategicRouteRealization:
    """One explicit route assessed against a selected target partition."""

    realization_id: str
    partition_id: str
    route_id: str
    route_tree_id: str
    route_tree: ReactionRouteTree
    route_cost: float
    reaction_count: int
    best_frontier: PartitionFrontierMatch
    frontier_states: tuple[LatentModuleState, ...]
    interface_realizations: tuple[InterfaceRealization, ...]
    dependency_edges: tuple[tuple[str, str], ...]
    status: str
    evidence_summary: RealizationEvidenceSummary
    verification_status: str
    warnings: tuple[str, ...]
    schema_version: str = PARTITION_REALIZATION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != PARTITION_REALIZATION_SCHEMA_VERSION:
            raise ValueError("unsupported strategic realization schema")
        if self.status not in PARTITION_REALIZATION_STATUSES - {"not_attempted"}:
            raise ValueError("unsupported strategic realization status")
        if self.route_tree.tree_id != self.route_tree_id:
            raise ValueError("route-tree identity mismatch")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible strategic realization."""

        return {
            "realization_id": self.realization_id,
            "partition_id": self.partition_id,
            "route_id": self.route_id,
            "route_tree_id": self.route_tree_id,
            "route_tree": self.route_tree.to_dict(),
            "route_cost": self.route_cost,
            "reaction_count": self.reaction_count,
            "best_frontier": self.best_frontier.to_dict(),
            "frontier_states": [item.to_dict() for item in self.frontier_states],
            "interface_realizations": [
                item.to_dict() for item in self.interface_realizations
            ],
            "dependency_edges": [list(item) for item in self.dependency_edges],
            "status": self.status,
            "evidence_summary": self.evidence_summary.to_dict(),
            "verification_status": self.verification_status,
            "warnings": list(self.warnings),
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "StrategicRouteRealization":
        """Reconstruct a strategic route realization."""

        return cls(
            realization_id=str(value["realization_id"]),
            partition_id=str(value["partition_id"]),
            route_id=str(value["route_id"]),
            route_tree_id=str(value["route_tree_id"]),
            route_tree=ReactionRouteTree.from_dict(dict(value["route_tree"])),
            route_cost=float(value["route_cost"]),
            reaction_count=int(value["reaction_count"]),
            best_frontier=PartitionFrontierMatch.from_dict(value["best_frontier"]),
            frontier_states=tuple(
                LatentModuleState.from_dict(item)
                for item in value.get("frontier_states") or ()
            ),
            interface_realizations=tuple(
                InterfaceRealization.from_dict(item)
                for item in value.get("interface_realizations") or ()
            ),
            dependency_edges=tuple(
                (str(item[0]), str(item[1]))
                for item in value.get("dependency_edges") or ()
            ),
            status=str(value["status"]),
            evidence_summary=RealizationEvidenceSummary.from_dict(
                value["evidence_summary"]
            ),
            verification_status=str(value["verification_status"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            schema_version=str(
                value.get("schema_version") or PARTITION_REALIZATION_SCHEMA_VERSION
            ),
        )


@dataclass(frozen=True)
class PartitionRealizationDiagnostics:
    """Bounded-search and projection counters for one request."""

    expanded_states: int
    generated_candidates: int
    rejected_cycles: int
    rejected_invalid_candidates: int
    dead_end_states: int
    stopped_by_expansion_limit: bool
    candidate_route_count: int
    projected_route_count: int
    projection_failure_count: int
    fully_realized_count: int
    partially_realized_count: int
    contradicted_count: int
    returned_realization_count: int

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible diagnostics."""

        return asdict(self)

    @classmethod
    def from_dict(
        cls,
        value: Mapping[str, Any],
    ) -> "PartitionRealizationDiagnostics":
        """Reconstruct realization diagnostics."""

        return cls(**{field: value[field] for field in cls.__dataclass_fields__})


@dataclass(frozen=True)
class PartitionRealizationResult:
    """Versioned result for one selected partition and bounded route search."""

    target_smiles: str
    partition: SyntheticPartition
    status: str
    realizations: tuple[StrategicRouteRealization, ...]
    config: PartitionRealizationConfig
    diagnostics: PartitionRealizationDiagnostics
    warnings: tuple[str, ...]
    policy_definition_id: str
    algorithm_version: str = PARTITION_REALIZATION_ALGORITHM_VERSION
    schema_version: str = PARTITION_REALIZATION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != PARTITION_REALIZATION_SCHEMA_VERSION:
            raise ValueError("unsupported partition realization result schema")
        if self.algorithm_version != PARTITION_REALIZATION_ALGORITHM_VERSION:
            raise ValueError("unsupported partition realization algorithm")
        if self.status not in PARTITION_REALIZATION_STATUSES - {"not_attempted"}:
            raise ValueError("unsupported partition realization result status")
        if self.partition.realization_status != self.status:
            raise ValueError("partition and result realization statuses disagree")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible realization result."""

        return {
            "target_smiles": self.target_smiles,
            "partition": self.partition.to_dict(),
            "status": self.status,
            "realizations": [item.to_dict() for item in self.realizations],
            "config": self.config.to_dict(),
            "diagnostics": self.diagnostics.to_dict(),
            "warnings": list(self.warnings),
            "policy_definition_id": self.policy_definition_id,
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionRealizationResult":
        """Reconstruct a serialized realization result."""

        return cls(
            target_smiles=str(value["target_smiles"]),
            partition=SyntheticPartition.from_dict(value["partition"]),
            status=str(value["status"]),
            realizations=tuple(
                StrategicRouteRealization.from_dict(item)
                for item in value.get("realizations") or ()
            ),
            config=PartitionRealizationConfig.from_dict(value["config"]),
            diagnostics=PartitionRealizationDiagnostics.from_dict(value["diagnostics"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            policy_definition_id=str(value["policy_definition_id"]),
            algorithm_version=str(value["algorithm_version"]),
            schema_version=str(value["schema_version"]),
        )


def validate_partition_realization_policy(value: Mapping[str, Any]) -> list[str]:
    """Return stable validation errors for the realization policy."""

    errors: list[str] = []
    if value.get("definition_id") != "synthetic_partition_realization_policy.v1":
        errors.append("invalid_definition_id")
    if value.get("schema_version") != "1.0":
        errors.append("invalid_schema_version")
    search = value.get("search")
    matching = value.get("matching")
    ranking = value.get("ranking")
    if not isinstance(search, Mapping):
        errors.append("missing_search")
        search = {}
    if not isinstance(matching, Mapping):
        errors.append("missing_matching")
        matching = {}
    if not isinstance(ranking, Mapping):
        errors.append("missing_ranking")
        ranking = {}
    try:
        PartitionRealizationConfig(
            max_depth=int(search.get("max_depth", 0)),
            molecular_weight_threshold=float(
                search.get("molecular_weight_threshold", 0.0)
            ),
            maximum_realizations=int(search.get("maximum_realizations", 0)),
            route_candidate_limit=int(search.get("route_candidate_limit", 0)),
            per_step_top_k=int(search.get("per_step_top_k", 0)),
            beam_width=int(search.get("beam_width", 0)),
            maximum_expansions=int(search.get("maximum_expansions", 0)),
            maximum_templates_to_apply=int(search.get("maximum_templates_to_apply", 0)),
            maximum_candidates_to_validate=int(
                search.get("maximum_candidates_to_validate", 0)
            ),
        )
    except (TypeError, ValueError):
        errors.append("invalid_search")
    try:
        threshold = float(matching.get("partial_interface_coverage_threshold", -1.0))
        if not 0.0 <= threshold <= 1.0:
            errors.append("invalid_partial_interface_coverage_threshold")
    except (TypeError, ValueError):
        errors.append("invalid_partial_interface_coverage_threshold")
    if tuple(ranking.get("status_order") or ()) != _STATUS_ORDER:
        errors.append("invalid_status_order")
    return errors


@lru_cache(maxsize=1)
def load_partition_realization_policy() -> PartitionRealizationPolicy:
    """Load and validate deterministic realization defaults."""

    value = json.loads(PARTITION_REALIZATION_POLICY_PATH.read_text("utf-8"))
    errors = validate_partition_realization_policy(value)
    if errors:
        raise ValueError("invalid partition realization policy: " + ", ".join(errors))
    search = value["search"]
    return PartitionRealizationPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        max_depth=int(search["max_depth"]),
        molecular_weight_threshold=float(search["molecular_weight_threshold"]),
        maximum_realizations=int(search["maximum_realizations"]),
        route_candidate_limit=int(search["route_candidate_limit"]),
        per_step_top_k=int(search["per_step_top_k"]),
        beam_width=int(search["beam_width"]),
        maximum_expansions=int(search["maximum_expansions"]),
        maximum_templates_to_apply=int(search["maximum_templates_to_apply"]),
        maximum_candidates_to_validate=int(search["maximum_candidates_to_validate"]),
        partial_interface_coverage_threshold=float(
            value["matching"]["partial_interface_coverage_threshold"]
        ),
        status_order=tuple(str(item) for item in value["ranking"]["status_order"]),
    )


def default_partition_realization_config() -> PartitionRealizationConfig:
    """Return a concrete configuration from the versioned policy."""

    policy = load_partition_realization_policy()
    return PartitionRealizationConfig(
        max_depth=policy.max_depth,
        molecular_weight_threshold=policy.molecular_weight_threshold,
        maximum_realizations=policy.maximum_realizations,
        route_candidate_limit=policy.route_candidate_limit,
        per_step_top_k=policy.per_step_top_k,
        beam_width=policy.beam_width,
        maximum_expansions=policy.maximum_expansions,
        maximum_templates_to_apply=policy.maximum_templates_to_apply,
        maximum_candidates_to_validate=policy.maximum_candidates_to_validate,
    )


def _module_sets(partition: SyntheticPartition) -> tuple[frozenset[int], ...]:
    return tuple(frozenset(module.target_atom_maps) for module in partition.modules)


def _partition_boundary_bonds(
    partition: SyntheticPartition,
) -> tuple[tuple[int, int, str], ...]:
    molecule = Chem.MolFromSmiles(partition.target_smiles)
    if molecule is None:
        raise ValueError("partition target could not be parsed")
    _, _, _, index_to_map = analyze_partition_target(partition.target_smiles)
    owner = {
        atom_map: module_index
        for module_index, module in enumerate(partition.modules)
        for atom_map in module.target_atom_maps
    }
    values = []
    for bond in molecule.GetBonds():
        left = index_to_map.get(bond.GetBeginAtomIdx())
        right = index_to_map.get(bond.GetEndAtomIdx())
        if left is None or right is None or owner[left] == owner[right]:
            continue
        values.append((min(left, right), max(left, right), str(bond.GetBondType())))
    return tuple(sorted(values))


def _frontier_boundary_bonds(
    frontier: RouteFrontierPartition,
) -> tuple[tuple[int, int, str], ...]:
    return tuple(
        sorted(
            {
                bond
                for interface in frontier.partition.interfaces
                for bond in interface.target_bonds
            }
        )
    )


def _atom_weighted_overlap(
    desired: Sequence[frozenset[int]],
    actual: Sequence[frozenset[int]],
) -> float:
    atom_count = sum(len(block) for block in desired)
    if atom_count == 0 or not actual:
        return 0.0
    desired_recall = (
        sum(
            max((len(block & other) for other in actual), default=0)
            for block in desired
        )
        / atom_count
    )
    actual_precision = (
        sum(
            max((len(block & other) for other in desired), default=0)
            for block in actual
        )
        / atom_count
    )
    if desired_recall + actual_precision == 0.0:
        return 0.0
    return round(
        2.0 * desired_recall * actual_precision / (desired_recall + actual_precision),
        8,
    )


def assess_partition_frontier_match(
    partition: SyntheticPartition,
    frontier: RouteFrontierPartition,
) -> PartitionFrontierMatch:
    """Compare one projected route frontier with a desired partition."""

    desired_modules = _module_sets(partition)
    actual_modules = _module_sets(frontier.partition)
    desired_set = frozenset(desired_modules)
    actual_set = frozenset(actual_modules)
    desired_bonds = set(_partition_boundary_bonds(partition))
    actual_bonds = set(_frontier_boundary_bonds(frontier))
    matched_bonds = desired_bonds & actual_bonds
    precision = (
        len(matched_bonds) / len(actual_bonds)
        if actual_bonds
        else 1.0
        if not desired_bonds
        else 0.0
    )
    recall = len(matched_bonds) / len(desired_bonds) if desired_bonds else 1.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return PartitionFrontierMatch(
        frontier_id=frontier.frontier_id,
        depth=frontier.depth,
        desired_k=partition.k,
        frontier_k=frontier.partition.k,
        exact_partition_match=desired_set == actual_set,
        matched_module_count=len(desired_set & actual_set),
        atom_weighted_module_overlap=_atom_weighted_overlap(
            desired_modules,
            actual_modules,
        ),
        boundary_precision=round(precision, 8),
        boundary_recall=round(recall, 8),
        boundary_f1=round(f1, 8),
        desired_boundary_count=len(desired_bonds),
        matched_boundary_count=len(matched_bonds),
        extra_boundary_count=len(actual_bonds - desired_bonds),
    )


def _match_sort_key(match: PartitionFrontierMatch) -> tuple[Any, ...]:
    return (
        not match.exact_partition_match,
        -match.boundary_f1,
        -match.boundary_recall,
        -match.atom_weighted_module_overlap,
        abs(match.frontier_k - match.desired_k),
        match.depth,
        match.frontier_id,
    )


def _best_frontier_match(
    partition: SyntheticPartition,
    projection: RoutePartitionProjection,
) -> tuple[PartitionFrontierMatch, RouteFrontierPartition]:
    if not projection.frontiers:
        raise ValueError("route projection has no complete frontiers")
    values = tuple(
        (assess_partition_frontier_match(partition, frontier), frontier)
        for frontier in projection.frontiers
    )
    return min(values, key=lambda item: _match_sort_key(item[0]))


def _tree_from_guidance_state(state: MultistepGuidanceState) -> ReactionRouteTree:
    root_id = next(
        (step.product_node_id for step in state.steps if step.depth == 1),
        digest("RMOL1", state.target_smiles, "root"),
    )
    return build_canonical_route_tree(
        state.target_smiles,
        root_id,
        state.steps,
        state.leaves,
    )


class PartitionSearchGuidance:
    """Symmetric route-state ordering toward one target partition."""

    definition_id = PARTITION_REALIZATION_ALGORITHM_VERSION

    def __init__(self, partition: SyntheticPartition) -> None:
        self.partition = partition
        self.desired_bonds = _partition_boundary_bonds(partition)
        self._projection_cache: dict[
            tuple[Any, ...], RoutePartitionProjection | None
        ] = {}

    @staticmethod
    def _state_key(state: MultistepGuidanceState) -> tuple[Any, ...]:
        return (
            tuple(step.candidate.proposed_reaction_smiles for step in state.steps),
            tuple(
                sorted(
                    (
                        leaf.route_node_id,
                        leaf.canonical_smiles,
                        leaf.depth,
                        leaf.terminal,
                        leaf.unresolved_reason or "",
                    )
                    for leaf in state.leaves
                )
            ),
        )

    def _project(
        self,
        state: MultistepGuidanceState,
    ) -> RoutePartitionProjection | None:
        key = self._state_key(state)
        if key not in self._projection_cache:
            try:
                self._projection_cache[key] = project_route_partitions(
                    _tree_from_guidance_state(state)
                )
            except ValueError:
                self._projection_cache[key] = None
        return self._projection_cache[key]

    def state_priority(
        self,
        state: MultistepGuidanceState,
    ) -> tuple[Any, ...]:
        projection = self._project(state)
        if projection is None or not projection.frontiers:
            return (1, 1, 1.0, 1.0, len(state.steps))
        match, _ = _best_frontier_match(self.partition, projection)
        return (
            0,
            not match.exact_partition_match,
            -match.boundary_recall,
            -match.atom_weighted_module_overlap,
            len(state.steps),
        )

    def select_leaf(
        self,
        state: MultistepGuidanceState,
        expandable_indices: tuple[int, ...],
    ) -> int:
        projection = self._project(state)
        ownership: dict[str, frozenset[int]] = {}
        if projection is not None and projection.frontiers:
            latest = max(projection.frontiers, key=lambda item: item.depth)
            ownership = {
                latent.source_occurrence_id: frozenset(latent.target_atom_maps)
                for latent in latest.latent_states
                if latent.source_occurrence_id is not None
            }
        desired_modules = _module_sets(self.partition)

        def priority(index: int) -> tuple[Any, ...]:
            leaf = state.leaves[index]
            atom_maps = ownership.get(leaf.route_node_id, frozenset())
            unresolved_bonds = sum(
                left in atom_maps and right in atom_maps
                for left, right, _ in self.desired_bonds
            )
            mixed_modules = sum(bool(atom_maps & block) for block in desired_modules)
            return (
                -unresolved_bonds,
                -mixed_modules,
                -(leaf.molecular_weight or 0.0),
                leaf.canonical_smiles,
                index,
            )

        return min(expandable_indices, key=priority)


def _dependency_edges(
    route: MultistepRetrosynthesisRoute,
) -> tuple[tuple[str, str], ...]:
    step_by_product = {step.product_node_id: step for step in route.steps}
    edges = {
        (child.step_id, parent.step_id)
        for parent in route.steps
        for precursor_id in parent.precursor_node_ids
        for child in (step_by_product.get(precursor_id),)
        if child is not None
    }
    return tuple(sorted(edges))


def _dependency_graph_valid(
    step_ids: Sequence[str],
    edges: Sequence[tuple[str, str]],
) -> bool:
    adjacency = {step_id: set() for step_id in step_ids}
    indegree = {step_id: 0 for step_id in step_ids}
    for source, target in edges:
        if source not in adjacency or target not in adjacency:
            return False
        if target in adjacency[source]:
            continue
        adjacency[source].add(target)
        indegree[target] += 1
    frontier = sorted(key for key, value in indegree.items() if value == 0)
    visited = 0
    while frontier:
        node = frontier.pop(0)
        visited += 1
        for target in sorted(adjacency[node]):
            indegree[target] -= 1
            if indegree[target] == 0:
                frontier.append(target)
                frontier.sort()
    return visited == len(indegree)


def _interface_realization(
    interface: StrategicInterface,
    route: MultistepRetrosynthesisRoute,
    projection: RoutePartitionProjection,
) -> InterfaceRealization:
    required = set(interface.target_bonds)
    matching_frontiers = tuple(
        frontier
        for frontier in projection.frontiers
        if required and required.issubset(set(_frontier_boundary_bonds(frontier)))
    )
    operator_steps = tuple(
        step
        for step in route.steps
        if (
            step.candidate.operator_id in interface.candidate_operator_ids
            or step.candidate.strategy_id in interface.candidate_strategy_ids
        )
    )
    if required:
        realized = bool(matching_frontiers)
        first_depth = min(
            (frontier.depth for frontier in matching_frontiers),
            default=None,
        )
        relevant_steps = operator_steps or tuple(
            step
            for step in route.steps
            if first_depth is not None and step.depth <= first_depth
        )
        kind = (
            "single_step_operator"
            if len(relevant_steps) == 1
            else "validated_route_path"
            if relevant_steps
            else "unresolved"
        )
    else:
        realized = bool(operator_steps)
        relevant_steps = operator_steps
        matching_frontiers = ()
        kind = "validated_unary_operator" if realized else "unresolved"
    matched = (
        required
        if realized
        else required.intersection(
            {
                bond
                for frontier in projection.frontiers
                for bond in _frontier_boundary_bonds(frontier)
            }
        )
    )
    warnings = () if realized else ("INTERFACE_NOT_REALIZED_IN_BOUNDED_ROUTE",)
    evidence_level = (
        "E3" if interface.evidence_level == "E4" else interface.evidence_level
    )
    return InterfaceRealization(
        interface_id=interface.interface_id,
        status="realized" if realized else "unresolved",
        realization_kind=kind,
        target_bonds=tuple(sorted(required)),
        matched_target_bonds=tuple(sorted(matched)),
        route_frontier_ids=tuple(
            frontier.frontier_id for frontier in matching_frontiers
        ),
        route_step_ids=tuple(sorted(step.step_id for step in relevant_steps)),
        operator_ids=tuple(
            sorted({step.candidate.operator_id for step in relevant_steps})
        ),
        strategy_ids=tuple(
            sorted(
                {
                    step.candidate.strategy_id
                    for step in relevant_steps
                    if step.candidate.strategy_id
                }
            )
        ),
        evidence_level=evidence_level,
        warnings=warnings,
    )


def _hard_verification_failed(
    report: RouteVerificationReport,
    *,
    route_solved: bool,
) -> bool:
    hard_gates = {
        "route_tree_integrity",
        "target_identity",
        "step_graph_consistency",
        "forward_validation",
    }
    if route_solved:
        hard_gates.add("chemistry_issue_screen")
    return any(
        gate.gate in hard_gates and gate.status == "fail" for gate in report.gates
    )


def _route_realization(
    partition: SyntheticPartition,
    route: MultistepRetrosynthesisRoute,
    projection: RoutePartitionProjection,
    policy: PartitionRealizationPolicy,
) -> StrategicRouteRealization:
    best_match, best_frontier = _best_frontier_match(partition, projection)
    interface_values = tuple(
        _interface_realization(interface, route, projection)
        for interface in partition.interfaces
    )
    realized_count = sum(item.status == "realized" for item in interface_values)
    interface_coverage = (
        realized_count / len(interface_values)
        if interface_values
        else best_match.boundary_recall
    )
    dependencies = _dependency_edges(route)
    dependency_valid = _dependency_graph_valid(
        tuple(step.step_id for step in route.steps),
        dependencies,
    )
    verification = verify_planned_route(route)
    contradicted = not dependency_valid or _hard_verification_failed(
        verification,
        route_solved=route.solved,
    )
    all_interfaces = realized_count == len(interface_values)
    if (
        not contradicted
        and route.solved
        and best_match.exact_partition_match
        and all_interfaces
    ):
        status = "fully_realized"
    elif contradicted:
        status = "contradicted"
    elif route.steps and (
        best_match.exact_partition_match
        or interface_coverage >= policy.partial_interface_coverage_threshold
        or best_match.boundary_recall >= policy.partial_interface_coverage_threshold
    ):
        status = "partially_realized"
    else:
        status = "unrealized_but_plausible"
    warnings = set(route.warnings) | set(projection.warnings)
    if not best_match.exact_partition_match:
        warnings.add("PARTITION_NOT_EXACTLY_EXPOSED")
    if not all_interfaces:
        warnings.add("PARTITION_INTERFACES_INCOMPLETE")
    if not route.solved:
        warnings.add("ROUTE_TERMINALS_INCOMPLETE")
    if partition.source_kind == "operator_combination_unrealized":
        warnings.add("SOURCE_PARTITION_WAS_UNREALIZED_HYPOTHESIS")
    if contradicted:
        warnings.add("HARD_ROUTE_CONTRADICTION")
    weakest = None
    if interface_values:
        weakest = max(
            (item.evidence_level for item in interface_values),
            key=_EVIDENCE_ORDER.index,
        )
    validated_steps = sum(
        step.candidate.forward_validation_status == "verified_signature"
        for step in route.steps
    )
    summary = RealizationEvidenceSummary(
        requested_interface_count=len(interface_values),
        realized_interface_count=realized_count,
        validated_interface_coverage=round(interface_coverage, 8),
        exact_partition_match=best_match.exact_partition_match,
        atom_weighted_module_overlap=best_match.atom_weighted_module_overlap,
        boundary_precision=best_match.boundary_precision,
        boundary_recall=best_match.boundary_recall,
        boundary_f1=best_match.boundary_f1,
        weakest_interface_evidence=weakest,
        validated_step_count=validated_steps,
        unsupported_step_count=len(route.steps) - validated_steps,
        dependency_graph_valid=dependency_valid,
        forward_order_exists=dependency_valid,
        terminal_leaf_count=sum(leaf.terminal for leaf in route.leaves),
        unresolved_leaf_count=sum(not leaf.terminal for leaf in route.leaves),
    )
    realization_id = digest(
        "SPREAL1",
        partition.partition_id,
        route.route_id,
        best_frontier.frontier_id,
        status,
        PARTITION_REALIZATION_SCHEMA_VERSION,
    )
    return StrategicRouteRealization(
        realization_id=realization_id,
        partition_id=partition.partition_id,
        route_id=route.route_id,
        route_tree_id=route.route_tree.tree_id,
        route_tree=route.route_tree,
        route_cost=route.route_cost,
        reaction_count=route.reaction_count,
        best_frontier=best_match,
        frontier_states=best_frontier.latent_states,
        interface_realizations=interface_values,
        dependency_edges=dependencies,
        status=status,
        evidence_summary=summary,
        verification_status=verification.status,
        warnings=tuple(sorted(warnings)),
    )


def _realization_sort_key(
    realization: StrategicRouteRealization,
    policy: PartitionRealizationPolicy,
) -> tuple[Any, ...]:
    return (
        policy.status_order.index(realization.status),
        not realization.best_frontier.exact_partition_match,
        -realization.evidence_summary.validated_interface_coverage,
        -realization.best_frontier.boundary_f1,
        -realization.best_frontier.atom_weighted_module_overlap,
        realization.route_cost,
        realization.reaction_count,
        realization.realization_id,
    )


def realize_synthetic_partition(
    partition: SyntheticPartition,
    library: GenericTemplateLibrary,
    literature_index: LiteratureLookup,
    *,
    config: PartitionRealizationConfig | None = None,
    precursor_realism_scorer: PrecursorRealismScorer | None = None,
    condition_evidence_evaluator: ConditionEvidenceEvaluator | None = None,
    route_action_policy: RouteActionPolicyModel | None = None,
    candidate_exclusions: tuple[RouteCandidateExclusion, ...] = (),
    expander: OneStepExpander | None = None,
) -> PartitionRealizationResult:
    """Search for validated routes that expose one selected partition."""

    report = validate_synthetic_partition(partition)
    if not report.valid:
        raise ValueError("invalid selected partition: " + ", ".join(report.issues))
    canonical, target_id, _, _ = analyze_partition_target(partition.target_smiles)
    if canonical != partition.target_smiles or target_id != partition.target_id:
        raise ValueError("selected partition target identity is inconsistent")
    active_config = config or default_partition_realization_config()
    policy = load_partition_realization_policy()
    guidance = PartitionSearchGuidance(partition)
    search = plan_multistep_routes(
        canonical,
        library,
        literature_index,
        max_depth=active_config.max_depth,
        molecular_weight_threshold=active_config.molecular_weight_threshold,
        top_k_routes=active_config.route_candidate_limit,
        per_step_top_k=active_config.per_step_top_k,
        beam_width=active_config.beam_width,
        max_expansions=active_config.maximum_expansions,
        max_templates_to_apply=active_config.maximum_templates_to_apply,
        max_candidates_to_validate=(active_config.maximum_candidates_to_validate),
        use_context=active_config.use_context,
        include_l0=active_config.include_l0,
        diversify=active_config.diversify,
        use_hierarchical_ranking=active_config.use_hierarchical_ranking,
        precursor_realism_scorer=precursor_realism_scorer,
        condition_evidence_evaluator=condition_evidence_evaluator,
        route_action_policy=route_action_policy,
        candidate_exclusions=candidate_exclusions,
        expander=expander,
        search_guidance=guidance,
    )
    route_by_id = {
        route.route_id: route for route in (*search.routes, *search.partial_routes)
    }
    realized_values = []
    projection_failures = 0
    for route in route_by_id.values():
        try:
            projection = project_route_partitions(route.route_tree)
            realized_values.append(
                _route_realization(partition, route, projection, policy)
            )
        except ValueError:
            projection_failures += 1
    ordered = tuple(
        sorted(
            realized_values,
            key=lambda item: _realization_sort_key(item, policy),
        )[: active_config.maximum_realizations]
    )
    if ordered:
        status = ordered[0].status
    else:
        status = "unrealized_but_plausible"
    warnings: set[str] = set()
    if not route_by_id:
        warnings.add("NO_BOUNDED_ROUTE_CANDIDATES")
    elif not ordered:
        warnings.add("NO_PROJECTABLE_ROUTE_CANDIDATES")
    if search.diagnostics.stopped_by_expansion_limit:
        warnings.add("REALIZATION_SEARCH_LIMIT_REACHED")
    if search.diagnostics.dead_end_states:
        warnings.add("OPERATOR_COVERAGE_INCOMPLETE")
    if projection_failures:
        warnings.add("ROUTE_ATOM_PROJECTION_FAILURES")
    if partition.source_kind == "operator_combination_unrealized":
        warnings.add("SOURCE_PARTITION_WAS_UNREALIZED_HYPOTHESIS")
    diagnostics = PartitionRealizationDiagnostics(
        expanded_states=search.diagnostics.expanded_states,
        generated_candidates=search.diagnostics.generated_candidates,
        rejected_cycles=search.diagnostics.rejected_cycles,
        rejected_invalid_candidates=(search.diagnostics.rejected_invalid_candidates),
        dead_end_states=search.diagnostics.dead_end_states,
        stopped_by_expansion_limit=(search.diagnostics.stopped_by_expansion_limit),
        candidate_route_count=len(route_by_id),
        projected_route_count=len(realized_values),
        projection_failure_count=projection_failures,
        fully_realized_count=sum(
            item.status == "fully_realized" for item in realized_values
        ),
        partially_realized_count=sum(
            item.status == "partially_realized" for item in realized_values
        ),
        contradicted_count=sum(
            item.status == "contradicted" for item in realized_values
        ),
        returned_realization_count=len(ordered),
    )
    return PartitionRealizationResult(
        target_smiles=canonical,
        partition=replace(partition, realization_status=status),
        status=status,
        realizations=ordered,
        config=active_config,
        diagnostics=diagnostics,
        warnings=tuple(sorted(warnings)),
        policy_definition_id=policy.definition_id,
    )


__all__ = [
    "PARTITION_REALIZATION_ALGORITHM_VERSION",
    "PARTITION_REALIZATION_POLICY_PATH",
    "PARTITION_REALIZATION_SCHEMA_VERSION",
    "InterfaceRealization",
    "PartitionFrontierMatch",
    "PartitionRealizationConfig",
    "PartitionRealizationDiagnostics",
    "PartitionRealizationPolicy",
    "PartitionRealizationResult",
    "PartitionSearchGuidance",
    "RealizationEvidenceSummary",
    "StrategicRouteRealization",
    "assess_partition_frontier_match",
    "default_partition_realization_config",
    "load_partition_realization_policy",
    "realize_synthetic_partition",
    "validate_partition_realization_policy",
]
