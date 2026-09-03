"""Auditable latent-state transitions and forward realization schedules."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Mapping

from .chemistry import digest
from .partition_projection import RoutePartitionProjection
from .route_contract import MoleculeOccurrenceNode, ReactionRouteTree
from .synthetic_partition import LatentModuleState, SyntheticPartition


LATENT_TRANSITION_SCHEMA_VERSION = "1.0"
LATENT_TRANSITION_ALGORITHM_VERSION = "route_latent_transition_graph.v1"


@dataclass(frozen=True)
class LatentStateTransition:
    """One validated physical step connecting target-derived latent states."""

    transition_id: str
    step_id: str
    retrosynthetic_depth: int
    input_state_ids: tuple[str, ...]
    output_state_id: str
    tactical_input_occurrence_ids: tuple[str, ...]
    transition_kind: str
    action_class: str
    operator_id: str | None
    evidence_status: str
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible transition."""

        return {
            "transition_id": self.transition_id,
            "step_id": self.step_id,
            "retrosynthetic_depth": self.retrosynthetic_depth,
            "input_state_ids": list(self.input_state_ids),
            "output_state_id": self.output_state_id,
            "tactical_input_occurrence_ids": list(
                self.tactical_input_occurrence_ids
            ),
            "transition_kind": self.transition_kind,
            "action_class": self.action_class,
            "operator_id": self.operator_id,
            "evidence_status": self.evidence_status,
            "warnings": list(self.warnings),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "LatentStateTransition":
        """Reconstruct a serialized transition."""

        return cls(
            transition_id=str(value["transition_id"]),
            step_id=str(value["step_id"]),
            retrosynthetic_depth=int(value["retrosynthetic_depth"]),
            input_state_ids=tuple(
                str(item) for item in value.get("input_state_ids") or ()
            ),
            output_state_id=str(value["output_state_id"]),
            tactical_input_occurrence_ids=tuple(
                str(item)
                for item in value.get("tactical_input_occurrence_ids") or ()
            ),
            transition_kind=str(value["transition_kind"]),
            action_class=str(value["action_class"]),
            operator_id=(
                str(value["operator_id"])
                if value.get("operator_id") is not None
                else None
            ),
            evidence_status=str(value["evidence_status"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class ForwardRealizationStage:
    """A deterministic layer of mutually independent forward transitions."""

    stage_index: int
    transition_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible schedule stage."""

        return {
            "stage_index": self.stage_index,
            "transition_ids": list(self.transition_ids),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ForwardRealizationStage":
        """Reconstruct a serialized schedule stage."""

        return cls(
            stage_index=int(value["stage_index"]),
            transition_ids=tuple(
                str(item) for item in value.get("transition_ids") or ()
            ),
        )


@dataclass(frozen=True)
class LatentRealizationGraph:
    """Target-provenance state graph plus one canonical forward schedule."""

    source_tree_id: str
    states: tuple[LatentModuleState, ...]
    transitions: tuple[LatentStateTransition, ...]
    dependency_edges: tuple[tuple[str, str], ...]
    forward_stages: tuple[ForwardRealizationStage, ...]
    validation_status: str
    warnings: tuple[str, ...] = ()
    algorithm_version: str = LATENT_TRANSITION_ALGORITHM_VERSION
    schema_version: str = LATENT_TRANSITION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != LATENT_TRANSITION_SCHEMA_VERSION:
            raise ValueError("unsupported latent-transition schema")
        if self.algorithm_version != LATENT_TRANSITION_ALGORITHM_VERSION:
            raise ValueError("unsupported latent-transition algorithm")
        if self.validation_status != "validated":
            raise ValueError("latent realization graph must be validated")
        state_ids = tuple(state.latent_state_id for state in self.states)
        transition_ids = tuple(
            transition.transition_id for transition in self.transitions
        )
        if len(state_ids) != len(set(state_ids)):
            raise ValueError("latent realization graph has duplicate states")
        if len(transition_ids) != len(set(transition_ids)):
            raise ValueError("latent realization graph has duplicate transitions")
        known_states = set(state_ids)
        state_by_id = {state.latent_state_id: state for state in self.states}
        for transition in self.transitions:
            if transition.output_state_id not in known_states or not set(
                transition.input_state_ids
            ).issubset(known_states):
                raise ValueError("latent transition references an unknown state")
            input_maps = [
                atom_map
                for state_id in transition.input_state_ids
                for atom_map in state_by_id[state_id].target_atom_maps
            ]
            output_maps = state_by_id[
                transition.output_state_id
            ].target_atom_maps
            if len(input_maps) != len(set(input_maps)):
                raise ValueError(
                    "latent precursor states duplicate target ownership"
                )
            if set(input_maps) != set(output_maps):
                raise ValueError("latent transition changes target-atom ownership")
        known_transitions = set(transition_ids)
        if len(self.dependency_edges) != len(set(self.dependency_edges)):
            raise ValueError("latent realization graph has duplicate dependencies")
        if any(
            source not in known_transitions or target not in known_transitions
            for source, target in self.dependency_edges
        ):
            raise ValueError("latent dependency references an unknown transition")
        producer_by_state: dict[str, str] = {}
        for transition in self.transitions:
            if transition.output_state_id in producer_by_state:
                raise ValueError("latent state has multiple producing transitions")
            producer_by_state[transition.output_state_id] = transition.transition_id
        expected_dependencies = {
            (producer_by_state[state_id], transition.transition_id)
            for transition in self.transitions
            for state_id in transition.input_state_ids
            if state_id in producer_by_state
        }
        if set(self.dependency_edges) != expected_dependencies:
            raise ValueError(
                "latent dependencies do not match state consumption"
            )
        scheduled = tuple(
            transition_id
            for stage in self.forward_stages
            for transition_id in stage.transition_ids
        )
        if len(scheduled) != len(set(scheduled)) or set(scheduled) != (
            known_transitions
        ):
            raise ValueError("forward schedule must cover every transition once")
        if tuple(stage.stage_index for stage in self.forward_stages) != tuple(
            range(1, len(self.forward_stages) + 1)
        ):
            raise ValueError("forward schedule stages must be contiguous")
        stage_by_transition = {
            transition_id: stage.stage_index
            for stage in self.forward_stages
            for transition_id in stage.transition_ids
        }
        if any(
            stage_by_transition[source] >= stage_by_transition[target]
            for source, target in self.dependency_edges
        ):
            raise ValueError("forward schedule contradicts a dependency")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible latent realization graph."""

        return {
            "source_tree_id": self.source_tree_id,
            "states": [state.to_dict() for state in self.states],
            "transitions": [item.to_dict() for item in self.transitions],
            "dependency_edges": [list(item) for item in self.dependency_edges],
            "forward_stages": [stage.to_dict() for stage in self.forward_stages],
            "validation_status": self.validation_status,
            "warnings": list(self.warnings),
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "LatentRealizationGraph":
        """Reconstruct and revalidate a serialized latent realization graph."""

        return cls(
            source_tree_id=str(value["source_tree_id"]),
            states=tuple(
                LatentModuleState.from_dict(item)
                for item in value.get("states") or ()
            ),
            transitions=tuple(
                LatentStateTransition.from_dict(item)
                for item in value.get("transitions") or ()
            ),
            dependency_edges=tuple(
                (str(item[0]), str(item[1]))
                for item in value.get("dependency_edges") or ()
            ),
            forward_stages=tuple(
                ForwardRealizationStage.from_dict(item)
                for item in value.get("forward_stages") or ()
            ),
            validation_status=str(value["validation_status"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            algorithm_version=str(
                value.get("algorithm_version")
                or LATENT_TRANSITION_ALGORITHM_VERSION
            ),
            schema_version=str(
                value.get("schema_version") or LATENT_TRANSITION_SCHEMA_VERSION
            ),
        )


def _route_nodes(
    root: MoleculeOccurrenceNode,
) -> tuple[MoleculeOccurrenceNode, ...]:
    values = [root]
    if root.reaction is not None:
        for child in root.reaction.children:
            values.extend(_route_nodes(child))
    return tuple(values)


def _selected_partition_state(
    state: LatentModuleState,
    partition: SyntheticPartition,
) -> LatentModuleState:
    target_maps = set(state.target_atom_maps)
    module_ids = tuple(
        module.module_id
        for module in partition.modules
        if set(module.target_atom_maps).issubset(target_maps)
    )
    annotations = list(state.state_annotations)
    if len(module_ids) > 1:
        annotations.append("JOINED_TARGET_MODULE_STATE")
    elif len(module_ids) == 1:
        annotations.append("SINGLE_TARGET_MODULE_STATE")
    else:
        annotations.append("PARTIAL_TARGET_MODULE_STATE")
    return replace(
        state,
        module_ids=module_ids,
        state_annotations=tuple(dict.fromkeys(annotations)),
    )


def _forward_stages(
    transition_ids: tuple[str, ...],
    edges: tuple[tuple[str, str], ...],
) -> tuple[ForwardRealizationStage, ...]:
    adjacency = {transition_id: set() for transition_id in transition_ids}
    indegree = {transition_id: 0 for transition_id in transition_ids}
    for source, target in edges:
        if target in adjacency[source]:
            continue
        adjacency[source].add(target)
        indegree[target] += 1
    stages = []
    frontier = tuple(sorted(key for key, value in indegree.items() if value == 0))
    visited = 0
    while frontier:
        stages.append(
            ForwardRealizationStage(
                stage_index=len(stages) + 1,
                transition_ids=frontier,
            )
        )
        next_frontier = []
        for source in frontier:
            visited += 1
            for target in sorted(adjacency[source]):
                indegree[target] -= 1
                if indegree[target] == 0:
                    next_frontier.append(target)
        frontier = tuple(sorted(next_frontier))
    if visited != len(transition_ids):
        raise ValueError("latent transition dependencies contain a cycle")
    return tuple(stages)


def build_latent_realization_graph(
    tree: ReactionRouteTree,
    projection: RoutePartitionProjection,
    partition: SyntheticPartition,
    *,
    action_classes: Mapping[str, str] | None = None,
) -> LatentRealizationGraph:
    """Build a validated state-transition graph from a projected route tree."""

    if projection.source_tree_id != tree.tree_id:
        raise ValueError("route projection does not belong to the supplied tree")
    if projection.unresolved_occurrence_ids:
        raise ValueError("cannot build transitions from unresolved atom projection")
    selected_states: dict[str, LatentModuleState] = {}
    for frontier in projection.frontiers:
        for state in frontier.latent_states:
            occurrence_id = state.source_occurrence_id
            if occurrence_id is None:
                continue
            selected_states[occurrence_id] = _selected_partition_state(
                state,
                partition,
            )
    nodes = _route_nodes(tree.root)
    transition_by_step: dict[str, LatentStateTransition] = {}
    transition_by_product: dict[str, str] = {}
    classes = action_classes or {}
    for node in nodes:
        reaction = node.reaction
        if reaction is None:
            continue
        output = selected_states.get(node.occurrence_id)
        if output is None:
            raise ValueError("route step product has no projected latent state")
        inputs = tuple(
            selected_states[child.occurrence_id]
            for child in reaction.children
            if child.occurrence_id in selected_states
        )
        if not inputs:
            raise ValueError("route step has no target-derived precursor state")
        input_maps = [atom_map for state in inputs for atom_map in state.target_atom_maps]
        if len(input_maps) != len(set(input_maps)):
            raise ValueError("latent precursor states duplicate target ownership")
        if set(input_maps) != set(output.target_atom_maps):
            raise ValueError("latent transition changes target-atom ownership")
        action_class = classes.get(reaction.step_id, "unclassified")
        transition_kind = (
            "ring_reorganization"
            if action_class == "ring_reorganization"
            else "assembly"
            if len(inputs) > 1
            else "state_change"
        )
        transition_id = digest(
            "LTRANS1",
            tree.tree_id,
            reaction.step_id,
            output.latent_state_id,
            *(sorted(state.latent_state_id for state in inputs)),
        )
        transition_by_step[reaction.step_id] = LatentStateTransition(
            transition_id=transition_id,
            step_id=reaction.step_id,
            retrosynthetic_depth=reaction.depth,
            input_state_ids=tuple(
                sorted(state.latent_state_id for state in inputs)
            ),
            output_state_id=output.latent_state_id,
            tactical_input_occurrence_ids=tuple(
                sorted(
                    child.occurrence_id
                    for child in reaction.children
                    if child.occurrence_id not in selected_states
                )
            ),
            transition_kind=transition_kind,
            action_class=action_class,
            operator_id=reaction.operator_id,
            evidence_status="validated_route_step_projection",
        )
        transition_by_product[node.occurrence_id] = transition_id
    dependency_edges = tuple(
        sorted(
            (
                transition_by_product[child.occurrence_id],
                transition_by_product[node.occurrence_id],
            )
            for node in nodes
            if node.reaction is not None
            for child in node.reaction.children
            if child.reaction is not None
        )
    )
    transitions = tuple(
        sorted(
            transition_by_step.values(),
            key=lambda item: (item.retrosynthetic_depth, item.step_id),
        )
    )
    transition_ids = tuple(item.transition_id for item in transitions)
    return LatentRealizationGraph(
        source_tree_id=tree.tree_id,
        states=tuple(
            sorted(
                selected_states.values(),
                key=lambda item: item.latent_state_id,
            )
        ),
        transitions=transitions,
        dependency_edges=dependency_edges,
        forward_stages=_forward_stages(transition_ids, dependency_edges),
        validation_status="validated",
    )


__all__ = [
    "LATENT_TRANSITION_ALGORITHM_VERSION",
    "LATENT_TRANSITION_SCHEMA_VERSION",
    "ForwardRealizationStage",
    "LatentRealizationGraph",
    "LatentStateTransition",
    "build_latent_realization_graph",
]
