"""Evidence-scoped connectivity observations and canonical shadow identities."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import defaultdict
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from .reaction_models import (
    AtomStateTransition,
    BondState,
    BondTransition,
    ConnectivityEditGraph,
    ConnectivityObservationScope,
    HydrogenDelta,
    ReactionAtomReference,
    ReactionEdit,
)


_LOCALIZED_BOND_UNITS = {
    "SINGLE": 1,
    "DOUBLE": 2,
    "TRIPLE": 3,
}
_CANONICALIZATION_PERMUTATION_LIMIT = 100_000

AtomProvenanceKey = Tuple[object, ...]


def atom_provenance_key(atom: ReactionAtomReference) -> AtomProvenanceKey:
    """Return a record-local key used to join references to the same atom."""
    if atom.atom_map_number is not None:
        return ("map", int(atom.atom_map_number))
    return (
        "position",
        atom.side,
        int(atom.component_index),
        int(atom.atom_index),
    )


def _canonical_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _digest(prefix: str, value: object, *, length: int = 32) -> str:
    payload = _canonical_json(value).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(payload).hexdigest()[:length]


def bond_state(order: Optional[str]) -> BondState:
    """Build a definite localized or typed bond/no-bond state."""
    return (
        BondState("no_bond", None)
        if order is None
        else BondState("bond", str(order).upper())
    )


def endpoint_absent_state() -> BondState:
    """Build a state for an endpoint absent from the reported product."""
    return BondState("endpoint_absent", None)


def unknown_bond_state() -> BondState:
    """Build a state whose connectivity is not established."""
    return BondState("unknown", None)


def _state_units(state: BondState) -> Optional[int]:
    if state.state_kind == "no_bond":
        return 0
    if state.state_kind != "bond":
        return None
    return _LOCALIZED_BOND_UNITS.get(str(state.order or "").upper())


def make_bond_transition(
    *,
    atom_1: ReactionAtomReference,
    atom_2: ReactionAtomReference,
    before_state: BondState,
    after_state: BondState,
    observation_scope: ConnectivityObservationScope,
    evidence: str,
    confidence: float,
) -> BondTransition:
    """Construct a transition and derive arithmetic only for localized bonds."""
    before_units = _state_units(before_state)
    after_units = _state_units(after_state)
    delta_units = (
        after_units - before_units
        if before_units is not None and after_units is not None
        else None
    )
    if atom_provenance_key(atom_2) < atom_provenance_key(atom_1):
        atom_1, atom_2 = atom_2, atom_1
    return BondTransition(
        atom_1=atom_1,
        atom_2=atom_2,
        before_state=before_state,
        after_state=after_state,
        delta_units=delta_units,
        observation_scope=observation_scope,
        evidence=evidence,
        confidence=confidence,
    )


def _atom_base_label(atom: ReactionAtomReference) -> Tuple[object, ...]:
    return (
        atom.element,
        int(atom.formal_charge),
        bool(atom.aromatic),
        atom.hybridization,
    )


def _bond_state_token(state: BondState) -> Tuple[str, str]:
    return (state.state_kind, state.order or "NONE")


def _transition_edge_label(transition: BondTransition) -> Tuple[object, ...]:
    return (
        _bond_state_token(transition.before_state),
        _bond_state_token(transition.after_state),
        (
            "NONE"
            if transition.delta_units is None
            else f"{transition.delta_units:+d}"
        ),
        transition.observation_scope,
    )


def _hydrogen_token(delta: HydrogenDelta) -> Tuple[object, ...]:
    return (
        int(delta.before_count),
        int(delta.after_count),
        int(delta.delta_count),
        delta.observation_scope,
    )


def _atom_state_token(transition: AtomStateTransition) -> Tuple[object, ...]:
    return (
        int(transition.before_formal_charge),
        int(transition.after_formal_charge),
        (
            "NONE"
            if transition.before_radical_electrons is None
            else str(transition.before_radical_electrons)
        ),
        (
            "NONE"
            if transition.after_radical_electrons is None
            else str(transition.after_radical_electrons)
        ),
        (
            "NONE"
            if transition.before_isotope is None
            else str(transition.before_isotope)
        ),
        (
            "NONE"
            if transition.after_isotope is None
            else str(transition.after_isotope)
        ),
        transition.observation_scope,
    )


def _assign_colors(payloads: Mapping[AtomProvenanceKey, object]) -> Dict[
    AtomProvenanceKey, int
]:
    unique = sorted({_canonical_json(value) for value in payloads.values()})
    color_by_token = {token: index for index, token in enumerate(unique)}
    return {
        node: color_by_token[_canonical_json(value)]
        for node, value in payloads.items()
    }


def _refine_colors(
    nodes: Sequence[AtomProvenanceKey],
    base_labels: Mapping[AtomProvenanceKey, object],
    annotations: Mapping[AtomProvenanceKey, object],
    adjacency: Mapping[
        AtomProvenanceKey,
        Sequence[Tuple[AtomProvenanceKey, Tuple[object, ...]]],
    ],
) -> Dict[AtomProvenanceKey, int]:
    payloads = {
        node: (base_labels[node], annotations.get(node, ())) for node in nodes
    }
    colors = _assign_colors(payloads)
    for _ in range(len(nodes) + 1):
        refined = {
            node: (
                base_labels[node],
                annotations.get(node, ()),
                colors[node],
                tuple(
                    sorted(
                        (
                            edge_label,
                            colors[neighbor],
                        )
                        for neighbor, edge_label in adjacency.get(node, ())
                    )
                ),
            )
            for node in nodes
        }
        new_colors = _assign_colors(refined)
        if new_colors == colors:
            return colors
        colors = new_colors
    return colors


def _component_encoding(
    node_order: Sequence[AtomProvenanceKey],
    base_labels: Mapping[AtomProvenanceKey, object],
    annotations: Mapping[AtomProvenanceKey, object],
    edges: Sequence[
        Tuple[AtomProvenanceKey, AtomProvenanceKey, Tuple[object, ...]]
    ],
) -> str:
    positions = {node: index for index, node in enumerate(node_order)}
    edge_tokens = tuple(
        sorted(
            (
                min(positions[left], positions[right]),
                max(positions[left], positions[right]),
                label,
            )
            for left, right, label in edges
        )
    )
    return _canonical_json(
        {
            "nodes": tuple(
                (base_labels[node], annotations.get(node, ()))
                for node in node_order
            ),
            "edges": edge_tokens,
        }
    )


def _fallback_component_encoding(
    nodes: Sequence[AtomProvenanceKey],
    colors: Mapping[AtomProvenanceKey, int],
    base_labels: Mapping[AtomProvenanceKey, object],
    annotations: Mapping[AtomProvenanceKey, object],
    edges: Sequence[
        Tuple[AtomProvenanceKey, AtomProvenanceKey, Tuple[object, ...]]
    ],
    adjacency: Mapping[
        AtomProvenanceKey,
        Sequence[Tuple[AtomProvenanceKey, Tuple[object, ...]]],
    ],
) -> str:
    node_tokens = tuple(
        sorted(
            (
                colors[node],
                base_labels[node],
                annotations.get(node, ()),
                tuple(
                    sorted(
                        (edge_label, colors[neighbor])
                        for neighbor, edge_label in adjacency.get(node, ())
                    )
                ),
            )
            for node in nodes
        )
    )
    edge_tokens = tuple(
        sorted(
            (
                min(colors[left], colors[right]),
                max(colors[left], colors[right]),
                label,
            )
            for left, right, label in edges
        )
    )
    return _canonical_json({"nodes": node_tokens, "edges": edge_tokens})


def _canonicalize_component(
    nodes: Sequence[AtomProvenanceKey],
    base_labels: Mapping[AtomProvenanceKey, object],
    annotations: Mapping[AtomProvenanceKey, object],
    edges: Sequence[
        Tuple[AtomProvenanceKey, AtomProvenanceKey, Tuple[object, ...]]
    ],
    adjacency: Mapping[
        AtomProvenanceKey,
        Sequence[Tuple[AtomProvenanceKey, Tuple[object, ...]]],
    ],
) -> Tuple[str, bool]:
    ordered_nodes = tuple(sorted(nodes))
    colors = _refine_colors(
        ordered_nodes,
        base_labels,
        annotations,
        adjacency,
    )
    cells: Dict[int, list[AtomProvenanceKey]] = defaultdict(list)
    for node in ordered_nodes:
        cells[colors[node]].append(node)
    ordered_cells = tuple(tuple(cells[color]) for color in sorted(cells))
    permutation_count = math.prod(math.factorial(len(cell)) for cell in ordered_cells)
    if permutation_count > _CANONICALIZATION_PERMUTATION_LIMIT:
        return (
            _fallback_component_encoding(
                ordered_nodes,
                colors,
                base_labels,
                annotations,
                edges,
                adjacency,
            ),
            True,
        )
    best: Optional[str] = None
    permutation_groups = tuple(
        tuple(itertools.permutations(cell)) for cell in ordered_cells
    )
    for selection in itertools.product(*permutation_groups):
        node_order = tuple(node for cell in selection for node in cell)
        encoded = _component_encoding(
            node_order,
            base_labels,
            annotations,
            edges,
        )
        if best is None or encoded < best:
            best = encoded
    if best is None:
        raise ValueError("cannot canonicalize an empty edit component")
    return best, False


def _transition_sort_key(transition: BondTransition) -> Tuple[object, ...]:
    return (
        atom_provenance_key(transition.atom_1),
        atom_provenance_key(transition.atom_2),
        _transition_edge_label(transition),
        transition.evidence,
    )


def _hydrogen_sort_key(delta: HydrogenDelta) -> Tuple[object, ...]:
    return (
        atom_provenance_key(delta.atom),
        _hydrogen_token(delta),
        delta.evidence,
    )


def _atom_state_sort_key(
    transition: AtomStateTransition,
) -> Tuple[object, ...]:
    return (
        atom_provenance_key(transition.reactant_atom),
        _atom_state_token(transition),
        transition.evidence,
    )


def build_connectivity_edit_graph(
    *,
    bond_transitions: Iterable[BondTransition],
    hydrogen_deltas: Iterable[HydrogenDelta] = (),
    atom_state_transitions: Iterable[AtomStateTransition] = (),
    evidence: str,
    confidence: float,
    warnings: Iterable[str] = (),
) -> ConnectivityEditGraph:
    """Build an order-invariant internal edit graph and shadow identity."""
    transitions = tuple(sorted(bond_transitions, key=_transition_sort_key))
    hydrogens = tuple(sorted(hydrogen_deltas, key=_hydrogen_sort_key))
    atom_states = tuple(
        sorted(atom_state_transitions, key=_atom_state_sort_key)
    )
    atoms: Dict[AtomProvenanceKey, ReactionAtomReference] = {}
    annotations: Dict[AtomProvenanceKey, list[Tuple[object, ...]]] = defaultdict(list)
    adjacency: Dict[
        AtomProvenanceKey,
        list[Tuple[AtomProvenanceKey, Tuple[object, ...]]],
    ] = defaultdict(list)
    edges: list[
        Tuple[AtomProvenanceKey, AtomProvenanceKey, Tuple[object, ...]]
    ] = []
    for transition in transitions:
        left = atom_provenance_key(transition.atom_1)
        right = atom_provenance_key(transition.atom_2)
        atoms.setdefault(left, transition.atom_1)
        atoms.setdefault(right, transition.atom_2)
        label = _transition_edge_label(transition)
        edges.append((left, right, label))
        adjacency[left].append((right, label))
        adjacency[right].append((left, label))
    for delta in hydrogens:
        node = atom_provenance_key(delta.atom)
        atoms.setdefault(node, delta.atom)
        annotations[node].append(("hydrogen", _hydrogen_token(delta)))
    for transition in atom_states:
        node = atom_provenance_key(transition.reactant_atom)
        atoms.setdefault(node, transition.reactant_atom)
        annotations[node].append(("atom_state", _atom_state_token(transition)))

    base_labels = {node: _atom_base_label(atom) for node, atom in atoms.items()}
    annotation_tokens = {
        node: tuple(sorted(values)) for node, values in annotations.items()
    }
    unvisited = set(atoms)
    components: list[Tuple[AtomProvenanceKey, ...]] = []
    while unvisited:
        seed = min(unvisited)
        stack = [seed]
        component: set[AtomProvenanceKey] = set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(
                neighbor
                for neighbor, _ in adjacency.get(node, ())
                if neighbor not in component
            )
        unvisited.difference_update(component)
        components.append(tuple(sorted(component)))

    component_encodings: list[str] = []
    overflow = False
    for component in components:
        component_set = set(component)
        component_edges = tuple(
            edge
            for edge in edges
            if edge[0] in component_set and edge[1] in component_set
        )
        encoding, component_overflow = _canonicalize_component(
            component,
            base_labels,
            annotation_tokens,
            component_edges,
            adjacency,
        )
        component_encodings.append(encoding)
        overflow = overflow or component_overflow
    component_encodings.sort()
    edit_component_keys = tuple(
        _digest("CEC1", json.loads(encoding)) for encoding in component_encodings
    )
    graph_payload = {
        "schema_version": "1.0",
        "components": tuple(json.loads(value) for value in component_encodings),
    }
    combined_warnings = set(warnings)
    if any(
        transition.delta_units is None
        and transition.before_state.state_kind in {"bond", "no_bond"}
        and transition.after_state.state_kind in {"bond", "no_bond"}
        and transition.before_state != transition.after_state
        for transition in transitions
    ):
        combined_warnings.add("UNSUPPORTED_BOND_DOMAIN")
    if overflow:
        combined_warnings.add("CONNECTIVITY_CANONICALIZATION_OVERFLOW")
    return ConnectivityEditGraph(
        bond_transitions=transitions,
        hydrogen_deltas=hydrogens,
        atom_state_transitions=atom_states,
        edit_component_keys=edit_component_keys,
        shadow_key=_digest("CEG1", graph_payload),
        evidence=evidence,
        confidence=confidence,
        warnings=tuple(sorted(combined_warnings)),
    )


def connectivity_graph_from_reaction_edits(
    edits: Iterable[ReactionEdit],
    *,
    observation_scope: ConnectivityObservationScope,
    evidence: str,
    confidence: float,
    hydrogen_before_counts: Optional[Mapping[AtomProvenanceKey, int]] = None,
    warnings: Iterable[str] = (),
) -> ConnectivityEditGraph:
    """Adapt current edits to the shadow model when original states are known."""
    hydrogen_before_counts = hydrogen_before_counts or {}
    transitions: list[BondTransition] = []
    hydrogen_changes: Dict[
        AtomProvenanceKey, Tuple[ReactionAtomReference, int, list[str]]
    ] = {}
    combined_warnings = set(warnings)
    for edit in edits:
        if edit.edit_type == "hydrogen_change":
            key = atom_provenance_key(edit.atom_1)
            direction = 1 if edit.new_order is not None else -1
            if key not in hydrogen_changes:
                hydrogen_changes[key] = (edit.atom_1, 0, [])
            atom, running, evidence_values = hydrogen_changes[key]
            evidence_values.append(edit.evidence)
            hydrogen_changes[key] = (atom, running + direction, evidence_values)
            continue
        if edit.atom_2 is None:
            combined_warnings.add("CONNECTIVITY_EDIT_ENDPOINT_MISSING")
            continue
        transitions.append(
            make_bond_transition(
                atom_1=edit.atom_1,
                atom_2=edit.atom_2,
                before_state=bond_state(edit.old_order),
                after_state=bond_state(edit.new_order),
                observation_scope=observation_scope,
                evidence=edit.evidence,
                confidence=edit.confidence,
            )
        )
    hydrogens: list[HydrogenDelta] = []
    for key, (atom, delta_count, evidence_values) in hydrogen_changes.items():
        before_count = hydrogen_before_counts.get(key)
        if before_count is None:
            combined_warnings.add("HYDROGEN_BEFORE_COUNT_UNAVAILABLE")
            continue
        after_count = before_count + delta_count
        if after_count < 0:
            combined_warnings.add("INVALID_PREDICTED_HYDROGEN_COUNT")
            continue
        hydrogens.append(
            HydrogenDelta(
                atom=atom,
                before_count=before_count,
                after_count=after_count,
                delta_count=delta_count,
                observation_scope=observation_scope,
                evidence=(
                    evidence_values[0]
                    if len(set(evidence_values)) == 1
                    else evidence
                ),
                confidence=confidence,
            )
        )
    return build_connectivity_edit_graph(
        bond_transitions=transitions,
        hydrogen_deltas=hydrogens,
        evidence=evidence,
        confidence=confidence,
        warnings=combined_warnings,
    )


def reaction_edits_from_connectivity_graph(
    graph: ConnectivityEditGraph,
) -> Tuple[ReactionEdit, ...]:
    """Project the stronger observation into the existing compatibility model."""
    edits: list[ReactionEdit] = []
    for transition in graph.bond_transitions:
        before = transition.before_state
        after = transition.after_state
        edit_type: Optional[str] = None
        old_order = before.order if before.state_kind == "bond" else None
        new_order = after.order if after.state_kind == "bond" else None
        if before.state_kind == "no_bond" and after.state_kind == "bond":
            edit_type = "formed"
        elif before.state_kind == "bond" and after.state_kind == "no_bond":
            edit_type = "broken"
        elif (
            before.state_kind == "bond"
            and after.state_kind == "bond"
            and before.order != after.order
        ):
            edit_type = "order_changed"
        elif (
            before.state_kind == "bond"
            and after.state_kind == "endpoint_absent"
        ):
            edit_type = "broken"
        if edit_type is None:
            continue
        edits.append(
            ReactionEdit(
                edit_type=edit_type,  # type: ignore[arg-type]
                atom_1=transition.atom_1,
                atom_2=transition.atom_2,
                old_order=old_order,
                new_order=new_order,
                evidence=transition.evidence,
                confidence=transition.confidence,
            )
        )
    for delta in graph.hydrogen_deltas:
        for _ in range(abs(delta.delta_count)):
            edits.append(
                ReactionEdit(
                    edit_type="hydrogen_change",
                    atom_1=delta.atom,
                    atom_2=None,
                    old_order="SINGLE" if delta.delta_count < 0 else None,
                    new_order="SINGLE" if delta.delta_count > 0 else None,
                    evidence=delta.evidence,
                    confidence=delta.confidence,
                )
            )
    return tuple(
        sorted(
            edits,
            key=lambda edit: (
                atom_provenance_key(edit.atom_1),
                (
                    atom_provenance_key(edit.atom_2)
                    if edit.atom_2 is not None
                    else ("hydrogen",)
                ),
                edit.edit_type,
                edit.old_order or "NONE",
                edit.new_order or "NONE",
            ),
        )
    )


__all__ = [
    "AtomProvenanceKey",
    "atom_provenance_key",
    "bond_state",
    "build_connectivity_edit_graph",
    "connectivity_graph_from_reaction_edits",
    "endpoint_absent_state",
    "make_bond_transition",
    "reaction_edits_from_connectivity_graph",
    "unknown_bond_state",
]
