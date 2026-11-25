"""Canonical before/after edit graphs for shared reaction-core projections.

Only observed center correspondence is joined across sides. Context atoms keep
their side identity; no correspondence is inferred for unmapped surroundings.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from typing import Any, Mapping

from rdkit import Chem

AtomPoint = tuple[int, int]


class SharedCoreProjectionLimit(ValueError):
    """The exact local graph exceeds a declared projection resource limit."""


@dataclass(frozen=True)
class SharedCoreAvailability:
    """Resource bounds on the incidence graph used by the exact canonicalizer."""

    definition_id: str
    maximum_center_distance_pairs: int
    maximum_incidence_vertices: int

    @property
    def fingerprint(self) -> str:
        """Bind derived coverage to its limits independently of chemistry identity."""
        return f"{self.definition_id}:{self.maximum_center_distance_pairs}:{self.maximum_incidence_vertices}"


@lru_cache(maxsize=1)
def shared_core_availability() -> SharedCoreAvailability:
    """Load the versioned local projection availability policy."""
    path = Path(__file__).with_name("definitions") / "shared_core_availability.v1.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.pop("schema_version") != "1.0" or data.get("definition_id") != "shared_core_availability.v1":
        raise ValueError("Invalid shared-core availability definition")
    for field in ("maximum_center_distance_pairs", "maximum_incidence_vertices"):
        if type(data.get(field)) is not int or data[field] < 1:
            raise ValueError(f"Invalid shared-core availability limit {field}")
    return SharedCoreAvailability(**data)


def atom_label(atom: Chem.Atom) -> list[Any]:
    """Retain chemical state and absolute stereochemistry independently of maps."""
    return [
        atom.GetSymbol(),
        atom.GetFormalCharge(),
        atom.GetIsAromatic(),
        str(atom.GetHybridization()),
        atom.GetIsotope(),
        atom.GetNumRadicalElectrons(),
        atom.GetTotalNumHs(),
        atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "",
    ]


def canonical_graph(labels: list[Any], edges: list[tuple[int, int, Any]]) -> str:
    """Canonicalize typed nodes/edges as a colored incidence graph."""
    encode = lambda value: json.dumps(value, sort_keys=True, separators=(",", ":"))
    encoded = [encode(["atom", label]) for label in labels]
    edge_labels = [encode(["edge", edge[2]]) for edge in edges]
    vocabulary = sorted(set(encoded + edge_labels))
    colors = {label: i + 1 for i, label in enumerate(vocabulary)}
    graph = Chem.RWMol()
    for label in encoded + edge_labels:
        atom = Chem.Atom(0)
        atom.SetIsotope(colors[label])
        graph.AddAtom(atom)
    for i, (left, right, _) in enumerate(edges, len(labels)):
        graph.AddBond(left, i, Chem.BondType.SINGLE)
        graph.AddBond(i, right, Chem.BondType.SINGLE)
    return encode([vocabulary, Chem.MolToSmiles(graph, canonical=True)])


def _point(state: Mapping[str, Any]) -> AtomPoint:
    return int(state["component_index"]), int(state["atom_index"])


def observed_edit_graph(
    reactants: Mapping[int, Chem.Mol],
    products: Mapping[int, Chem.Mol],
    transitions: Mapping[AtomPoint, Mapping[str, Any]],
    edits: tuple[Mapping[str, Any], ...],
    *,
    radius: int,
    limit: int,
    ports: Mapping[AtomPoint, list[str]],
    removed: frozenset[AtomPoint],
    generalized: frozenset[int],
    product_only: bool = False,
    actual_ports: bool = False,
) -> str:
    """Project protected edits, states, functional units and center topology.

    Qualified ports may remove only previously validated realization edits and
    departing atoms. Every other edit and before/after atom state stays explicit.
    Pairwise molecular distances retain ring-closure size and multi-event geometry
    even when local shells no longer connect the centers visually.
    """
    labels: list[Any] = []
    edges: list[tuple[int, int, Any]] = []
    centers: dict[AtomPoint, int] = {}
    side_centers: dict[str, dict[AtomPoint, int]] = {"before": {}, "after": {}}
    for point, transition in sorted(transitions.items()):
        if point in removed:
            continue
        before, after = transition.get("before_state"), transition.get("after_state")
        if product_only and not after:
            continue
        state_labels = []
        for side, state, molecules in (
            ("before", before, reactants),
            ("after", after, products),
        ):
            if product_only and side == "before":
                continue
            label = None
            if state:
                component, index = _point(state)
                atom = molecules[component].GetAtomWithIdx(index)
                label = atom_label(atom)
                if side == "before" and point in ports:
                    label[6] = "qualified_port"
                # Membership in a ring is chemical context; exact scaffold
                # membership is not a required routing key.
                label.append(atom.IsInRing())
                side_centers[side][_point(state)] = len(labels)
            state_labels.append(label)
        centers[point] = len(labels)
        labels.append(["center", *state_labels])

    availability = shared_core_availability()
    center_pairs = sum(
        count * (count - 1) // 2
        for points in side_centers.values()
        for count in Counter(component for component, _ in points).values()
    )
    if center_pairs > availability.maximum_center_distance_pairs:
        raise SharedCoreProjectionLimit("LOCAL_CORE_CENTER_DISTANCE_PAIR_LIMIT")
    for side, molecules in (("before", reactants), ("after", products)):
        if product_only and side == "before":
            continue
        for component, molecule in sorted(molecules.items()):
            anchors = {
                i: node for (c, i), node in side_centers[side].items() if c == component
            }
            if not anchors:
                continue
            excluded = (
                {i for c, i in removed if c == component} if side == "before" else set()
            )
            selected = set(anchors)
            for _ in range(radius):
                selected |= {
                    n.GetIdx()
                    for i in selected
                    for n in molecule.GetAtomWithIdx(i).GetNeighbors()
                } - excluded
            while True:
                expanded = (
                    selected
                    | {
                        bond.GetOtherAtomIdx(i)
                        for i in selected
                        for bond in molecule.GetAtomWithIdx(i).GetBonds()
                        if bond.GetBondType()
                        in {Chem.BondType.DOUBLE, Chem.BondType.TRIPLE}
                    }
                    - excluded
                )
                if len(expanded) > limit:
                    raise SharedCoreProjectionLimit("CORE_SIZE_LIMIT")
                if expanded == selected:
                    break
                selected = expanded
            nodes = dict(anchors)
            for i in sorted(selected - set(anchors)):
                nodes[i] = len(labels)
                labels.append([side, atom_label(molecule.GetAtomWithIdx(i))])
            for bond in molecule.GetBonds():
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if a in nodes and b in nodes:
                    edges.append(
                        (
                            nodes[a],
                            nodes[b],
                            [side, str(bond.GetBondType()), str(bond.GetStereo())],
                        )
                    )
            for a, b in combinations(sorted(anchors), 2):
                distance = len(Chem.GetShortestPath(molecule, a, b)) - 1
                if distance < 0:
                    raise ValueError("DISCONNECTED_COMPONENT")
                edges.append(
                    (anchors[a], anchors[b], [side, "center_distance", distance])
                )

    if not product_only:
        for i, edit in enumerate(edits):
            if i in generalized:
                continue
            references = [
                ref for ref in (edit.get("atom_1"), edit.get("atom_2")) if ref
            ]
            node = len(labels)
            labels.append(
                [
                    "edit",
                    edit["edit_type"],
                    edit.get("old_order"),
                    edit.get("new_order"),
                ]
            )
            for ref in references:
                edges.append((node, centers[_point(ref)], "endpoint"))
        for point in sorted(ports):
            node = len(labels)
            labels.append(
                ["observed_realization_port", sorted(ports[point])]
                if actual_ports
                else "qualified_realization_port"
            )
            edges.append((node, centers[point], "attachment"))
    # Every labeled edge becomes a separate vertex in the incidence graph.
    if len(labels) + len(edges) > availability.maximum_incidence_vertices:
        raise SharedCoreProjectionLimit("LOCAL_CORE_INCIDENCE_VERTEX_LIMIT")
    return canonical_graph(labels, edges)
