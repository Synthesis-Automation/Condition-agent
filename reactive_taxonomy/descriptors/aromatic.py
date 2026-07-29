"""Aromatic ring-system identity, positions, sterics, and electronics."""

from __future__ import annotations

from collections import Counter, deque
from typing import Any, Iterable, Tuple

from .common import branch_contribution, class_from_upper_bins, shortest_distance
from .models import (
    AromaticContextDescriptor,
    AromaticHeteroatom,
    ElectronicContribution,
    StericContribution,
)
from .registry import aromatic_system_rules, descriptor_rules


def _aromatic_rings(mol: Any) -> Tuple[Tuple[int, ...], ...]:
    return tuple(
        tuple(int(value) for value in ring)
        for ring in mol.GetRingInfo().AtomRings()
        if all(mol.GetAtomWithIdx(int(value)).GetIsAromatic() for value in ring)
    )


def aromatic_ring_system(
    mol: Any, center: int
) -> Tuple[Tuple[int, ...], Tuple[Tuple[int, ...], ...]]:
    """Return the fused aromatic system containing the center atom."""
    rings = _aromatic_rings(mol)
    seeds = {index for index, ring in enumerate(rings) if center in ring}
    if not seeds:
        return (), ()
    selected = set(seeds)
    queue = deque(seeds)
    while queue:
        current = queue.popleft()
        current_atoms = set(rings[current])
        for index, ring in enumerate(rings):
            if index in selected:
                continue
            if len(current_atoms & set(ring)) >= 2:
                selected.add(index)
                queue.append(index)
    system_rings = tuple(rings[index] for index in sorted(selected))
    system_atoms = tuple(sorted({atom for ring in system_rings for atom in ring}))
    return system_atoms, system_rings


def aromatic_atom_role(atom: Any) -> str:
    """Classify an aromatic heteroatom by lone-pair participation."""
    if atom.GetSymbol() == "N":
        if (
            atom.GetTotalNumHs() > 0
            or atom.GetDegree() >= 3
            or atom.GetFormalCharge() < 0
        ):
            return "pyrrole_like"
        return "pyridine_like"
    if atom.GetSymbol() in {"O", "S"}:
        return "pyrrole_like"
    return "other"


def _cyclic_distance(ring: Tuple[int, ...], atom_1: int, atom_2: int) -> int:
    first = ring.index(atom_1)
    second = ring.index(atom_2)
    direct = abs(first - second)
    return min(direct, len(ring) - direct)


def _heteroatom_records(
    mol: Any,
    center: int,
    system_atoms: Tuple[int, ...],
    system_rings: Tuple[Tuple[int, ...], ...],
) -> Tuple[AromaticHeteroatom, ...]:
    anchor_rings = tuple(ring for ring in system_rings if center in ring)
    records = []
    for atom_index in system_atoms:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() == 6:
            continue
        same_rings = tuple(
            ring
            for ring in anchor_rings
            if atom_index in ring
        )
        if atom_index == center:
            distance = 0
            relation = "anchor"
            same_anchor_ring = True
        elif same_rings:
            distance = min(
                _cyclic_distance(ring, center, atom_index) for ring in same_rings
            )
            max_ring_size = max(len(ring) for ring in same_rings)
            if distance == 1:
                relation = "ortho"
            elif distance == 2:
                relation = "meta"
            elif distance == 3 and max_ring_size == 6:
                relation = "para"
            else:
                relation = "remote"
            same_anchor_ring = True
        else:
            distance = shortest_distance(mol, center, (atom_index,)) or 0
            relation = "fused_other_ring"
            same_anchor_ring = False
        records.append(
            AromaticHeteroatom(
                atom_index=atom_index,
                element=atom.GetSymbol(),
                formal_charge=int(atom.GetFormalCharge()),
                aromatic_role=aromatic_atom_role(atom),  # type: ignore[arg-type]
                ring_distance_from_anchor=distance,
                positional_relation=relation,  # type: ignore[arg-type]
                same_anchor_ring=same_anchor_ring,
            )
        )
    return tuple(
        sorted(
            records,
            key=lambda item: (
                item.ring_distance_from_anchor,
                item.element,
                item.aromatic_role,
                item.formal_charge,
                item.atom_index,
            ),
        )
    )


def _composition(mol: Any, atoms: Iterable[int]) -> str:
    counts = Counter(mol.GetAtomWithIdx(index).GetSymbol() for index in atoms)
    return "".join(
        f"{element}{counts[element]}"
        for element in sorted(counts, key=lambda value: (value != "C", value))
    )


def _ring_family(
    mol: Any,
    system_atoms: Tuple[int, ...],
    system_rings: Tuple[Tuple[int, ...], ...],
    heteroatoms: Tuple[AromaticHeteroatom, ...],
) -> str:
    rules = aromatic_system_rules()
    if len(system_rings) == 1:
        ring = system_rings[0]
        composition = _composition(mol, ring)
        key = f"{len(ring)}:{composition}"
        if (
            len(heteroatoms) == 1
            and len(ring) == 5
            and heteroatoms[0].element == "N"
        ):
            key += f":{heteroatoms[0].aromatic_role}"
        elif len(heteroatoms) == 2 and len(ring) in {5, 6}:
            distance = _cyclic_distance(
                ring, heteroatoms[0].atom_index, heteroatoms[1].atom_index
            )
            relation = (
                "adjacent"
                if distance == 1
                else "opposite"
                if len(ring) == 6 and distance == 3
                else "separated"
            )
            key += f":{relation}"
        value = (rules.get("monocyclic_families") or {}).get(key)
        if value:
            return str(value)
    composition = _composition(mol, system_atoms)
    base = f"rings{len(system_rings)}:atoms{len(system_atoms)}:{composition}"
    roles = {item.aromatic_role for item in heteroatoms}
    candidates = [base]
    if len(roles) == 1:
        candidates.insert(0, f"{base}:{next(iter(roles))}")
    fused = rules.get("fused_families") or {}
    for key in candidates:
        if key in fused:
            return str(fused[key])
    fallback = rules.get("fallbacks") or {}
    if not heteroatoms:
        return str(fallback.get("carbocyclic") or "other_carbocyclic_aromatic")
    if len(system_rings) > 1:
        return str(fallback.get("mixed_fused") or "ambiguous_fused_system")
    return str(fallback.get("heteroaromatic") or "other_heteroaromatic")


def _ortho_sterics(
    mol: Any,
    center: int,
    system_atoms: Tuple[int, ...],
    system_rings: Tuple[Tuple[int, ...], ...],
) -> Tuple[int, int, str, float, Tuple[StericContribution, ...]]:
    anchor_ring_atoms = {
        atom
        for ring in system_rings
        if center in ring
        for atom in ring
    }
    ortho_atoms = tuple(
        sorted(
            neighbor.GetIdx()
            for neighbor in mol.GetAtomWithIdx(center).GetNeighbors()
            if neighbor.GetIdx() in anchor_ring_atoms and neighbor.GetIsAromatic()
        )
    )
    system_set = set(system_atoms)
    contributions = []
    occupied = 0
    for ortho in ortho_atoms:
        atom = mol.GetAtomWithIdx(ortho)
        external = tuple(
            neighbor.GetIdx()
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() > 1
            and neighbor.GetIdx() not in system_set
            and neighbor.GetIdx() != center
        )
        if external:
            occupied += 1
        for origin in external:
            contributions.append(
                branch_contribution(
                    mol,
                    origin,
                    blocked=system_set,
                    relation="ortho_substituent",
                    radius=3,
                )
            )
    raw = sum(max(0.35, value.score) for value in contributions)
    score = round(min(1.0, raw / max(1.0, len(ortho_atoms) * 0.5)), 3)
    bins = descriptor_rules()["steric"]["ortho_burden_bins"]
    burden = class_from_upper_bins(score, bins)
    return occupied, len(ortho_atoms), burden, score, tuple(
        sorted(
            contributions,
            key=lambda item: (item.relation, item.origin_atom_index),
        )
    )


def build_aromatic_context(
    mol: Any, center: int
) -> Tuple[
    AromaticContextDescriptor,
    Tuple[StericContribution, ...],
    Tuple[ElectronicContribution, ...],
]:
    """Build aromatic topology plus steric and intrinsic electronic evidence."""
    system_atoms, system_rings = aromatic_ring_system(mol, center)
    if not system_atoms:
        raise ValueError("aromatic context requires an aromatic anchor ring")
    heteroatoms = _heteroatom_records(
        mol, center, system_atoms, system_rings
    )
    ring_has_hetero = [
        any(mol.GetAtomWithIdx(index).GetAtomicNum() != 6 for index in ring)
        for ring in system_rings
    ]
    if any(ring_has_hetero) and not all(ring_has_hetero):
        system_class = "mixed_fused"
    elif any(ring_has_hetero):
        system_class = "heteroaromatic"
    else:
        system_class = "carbocyclic"
    occupied, capacity, burden, burden_score, steric = _ortho_sterics(
        mol, center, system_atoms, system_rings
    )
    weights = descriptor_rules()["electronic"]["aromatic_heteroatom_weights"]
    electronic = []
    for record in heteroatoms:
        if record.atom_index == center:
            continue
        key = f"{record.element}:{record.aromatic_role}"
        base = float(weights.get(key, weights.get("other", 0.0)))
        value = max(-1.0, min(1.0, base / max(1, record.ring_distance_from_anchor)))
        if value == 0.0:
            continue
        electronic.append(
            ElectronicContribution(
                source_id=f"aromatic_heteroatom:{record.element}",
                effect="withdrawing" if value > 0 else "donating",
                pathway="aromatic_intrinsic",
                positional_relation=record.positional_relation,
                contribution=round(value, 3),
                atom_indices=(record.atom_index,),
            )
        )
    context = AromaticContextDescriptor(
        context_kind="aromatic",
        system_class=system_class,  # type: ignore[arg-type]
        ring_family=_ring_family(
            mol, system_atoms, system_rings, heteroatoms
        ),
        ring_sizes=tuple(sorted(len(ring) for ring in system_rings)),
        aromatic_ring_count=len(system_rings),
        fused=len(system_rings) > 1,
        anchor_in_ring=True,
        heteroatoms=heteroatoms,
        ortho_occupancy_count=occupied,
        ortho_capacity=capacity,
        ortho_burden_class=burden,  # type: ignore[arg-type]
        ortho_burden_score=burden_score,
    )
    return (
        context,
        steric,
        tuple(
            sorted(
                electronic,
                key=lambda item: (
                    item.positional_relation,
                    item.source_id,
                    item.atom_indices,
                ),
            )
        ),
    )


__all__ = [
    "aromatic_atom_role",
    "aromatic_ring_system",
    "build_aromatic_context",
]
