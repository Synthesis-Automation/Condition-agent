"""Shared deterministic graph utilities for reactivity profiles."""

from __future__ import annotations

from collections import deque
from typing import Any, Iterable, Tuple

from .models import DescriptorEvidence, ReactiveCenterProfile, StericContribution


def atom_ring_sizes(mol: Any, atom_index: int) -> Tuple[int, ...]:
    """Return sorted ring sizes containing an atom."""
    return tuple(
        sorted(
            {
                len(ring)
                for ring in mol.GetRingInfo().AtomRings()
                if int(atom_index) in ring
            }
        )
    )


def shortest_distance(mol: Any, start: int, targets: Iterable[int]) -> int | None:
    """Return the shortest graph distance from one atom to target atoms."""
    from rdkit import Chem

    values = []
    for target in targets:
        if int(start) == int(target):
            values.append(0)
            continue
        try:
            values.append(len(Chem.GetShortestPath(mol, int(start), int(target))) - 1)
        except Exception:
            continue
    return min(values) if values else None


def bounded_branch_atoms(
    mol: Any,
    origin: int,
    blocked: Iterable[int],
    radius: int,
) -> Tuple[int, ...]:
    """Return heavy atoms in one bounded branch without crossing blocked atoms."""
    blocked_set = {int(value) for value in blocked}
    visited = set(blocked_set)
    visited.discard(int(origin))
    queue = deque([(int(origin), 0)])
    atoms = []
    while queue:
        atom_index, distance = queue.popleft()
        if atom_index in visited or distance > radius:
            continue
        visited.add(atom_index)
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() > 1:
            atoms.append(atom_index)
        if distance == radius:
            continue
        for neighbor in atom.GetNeighbors():
            index = neighbor.GetIdx()
            if index not in visited and neighbor.GetAtomicNum() > 1:
                queue.append((index, distance + 1))
    return tuple(sorted(atoms))


def branch_contribution(
    mol: Any,
    origin: int,
    blocked: Iterable[int],
    *,
    relation: str,
    radius: int = 3,
) -> StericContribution:
    """Build a normalized graph burden for one approach-side branch."""
    atoms = bounded_branch_atoms(mol, origin, blocked, radius)
    atom_set = set(atoms)
    branch_count = sum(
        max(
            0,
            sum(
                neighbor.GetAtomicNum() > 1
                and neighbor.GetIdx() in atom_set
                for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
            )
            - 2,
        )
        for index in atoms
    )
    score = min(1.0, (len(atoms) + 0.75 * branch_count) / 8.0)
    return StericContribution(
        origin_atom_index=int(origin),
        relation=relation,
        heavy_atom_count=len(atoms),
        branch_count=branch_count,
        score=round(score, 3),
    )


def carbon_substitution(atom: Any, *, excluded_atom: int | None = None) -> str:
    """Classify a carbon by attached carbon count."""
    carbon_neighbors = sum(
        neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != excluded_atom
        for neighbor in atom.GetNeighbors()
    )
    return {0: "methyl", 1: "primary", 2: "secondary"}.get(
        carbon_neighbors, "tertiary"
    )


def heteroatom_substitution(atom: Any) -> str:
    """Classify N/O/S/P substitution from heavy-atom attachments."""
    heavy = sum(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors())
    element = atom.GetSymbol()
    if element == "N":
        return {
            0: "ammonia",
            1: "primary",
            2: "secondary",
            3: "tertiary",
            4: "quaternary",
        }.get(heavy, "hypercoordinate")
    return f"{heavy}_coordinate"


def conjugation_class(atom: Any) -> str:
    """Return a small graph-derived conjugation vocabulary."""
    if atom.GetIsAromatic():
        return "aromatic"
    if any(neighbor.GetIsAromatic() for neighbor in atom.GetNeighbors()):
        return "aryl_conjugated"
    if any(
        str(bond.GetBondType()).upper() in {"DOUBLE", "TRIPLE", "AROMATIC"}
        for bond in atom.GetBonds()
    ):
        return "pi_conjugated"
    for neighbor in atom.GetNeighbors():
        if any(
            str(bond.GetBondType()).upper() in {"DOUBLE", "TRIPLE", "AROMATIC"}
            for bond in neighbor.GetBonds()
            if bond.GetOtherAtomIdx(neighbor.GetIdx()) != atom.GetIdx()
        ):
            return "adjacent_pi_system"
    return "saturated"


def reactive_center_profile(
    mol: Any,
    center: int,
    *,
    lone_pair_class: str | None,
    lone_pair_availability: str | None,
    acidity_class: str | None,
) -> ReactiveCenterProfile:
    """Build the context-independent atom-state profile."""
    atom = mol.GetAtomWithIdx(center)
    if atom.GetAtomicNum() == 6 and str(atom.GetHybridization()) in {"SP3", "S"}:
        substitution = carbon_substitution(atom)
    elif atom.GetAtomicNum() in {7, 8, 15, 16}:
        substitution = heteroatom_substitution(atom)
    else:
        substitution = None
    return ReactiveCenterProfile(
        atom_index=center,
        element=atom.GetSymbol(),
        hybridization=str(atom.GetHybridization()),
        formal_charge=int(atom.GetFormalCharge()),
        radical_electrons=int(atom.GetNumRadicalElectrons()),
        hydrogen_count=int(atom.GetTotalNumHs()),
        heavy_atom_attachment_count=sum(
            neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors()
        ),
        aromatic=bool(atom.GetIsAromatic()),
        in_ring=bool(atom.IsInRing()),
        ring_sizes=atom_ring_sizes(mol, center),
        substitution_class=substitution,
        conjugation_class=conjugation_class(atom),
        lone_pair_class=lone_pair_class,
        lone_pair_availability=lone_pair_availability,
        acidity_class=acidity_class,
        evidence=DescriptorEvidence(
            source="molecular_graph",
            method="reactive_center_graph_v1",
            confidence=1.0,
            contributing_atom_indices=(center,),
        ),
    )


def class_from_upper_bins(value: float, bins: dict[str, Any]) -> str:
    """Return the first ordered class whose upper bound contains a value."""
    for name, upper in bins.items():
        if value <= float(upper):
            return str(name)
    return str(next(reversed(bins)))


def electronic_class(value: float, bins: dict[str, Any]) -> str:
    """Return a signed electronic class from ordered boundary definitions."""
    ordered = tuple((str(name), float(boundary)) for name, boundary in bins.items())
    for name, boundary in ordered:
        if value <= boundary:
            return name
    return ordered[-1][0]


__all__ = [
    "atom_ring_sizes",
    "bounded_branch_atoms",
    "branch_contribution",
    "carbon_substitution",
    "class_from_upper_bins",
    "conjugation_class",
    "electronic_class",
    "heteroatom_substitution",
    "reactive_center_profile",
    "shortest_distance",
]
