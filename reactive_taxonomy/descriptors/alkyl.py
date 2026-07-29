"""Deterministic alkyl-center context descriptors."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from .common import branch_contribution, carbon_substitution
from .models import (
    AlkylContextDescriptor,
    ElectronicContribution,
    StericContribution,
)


def _adjacent_pi(atom: Any, order: str) -> bool:
    return any(
        str(bond.GetBondType()).upper() == order
        for neighbor in atom.GetNeighbors()
        for bond in neighbor.GetBonds()
        if bond.GetOtherAtomIdx(neighbor.GetIdx()) != atom.GetIdx()
    )


def build_alkyl_context(
    mol: Any,
    center: int,
    *,
    excluded_atoms: Iterable[int] = (),
) -> Tuple[
    AlkylContextDescriptor,
    Tuple[StericContribution, ...],
    Tuple[ElectronicContribution, ...],
]:
    """Build substitution, branching, activation, and accessibility evidence."""
    atom = mol.GetAtomWithIdx(center)
    if atom.GetAtomicNum() != 6:
        raise ValueError("alkyl context requires a carbon center")
    excluded = {int(value) for value in excluded_atoms}
    carbon_neighbors = tuple(
        neighbor
        for neighbor in atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in excluded
    )
    beta_branch_count = sum(
        max(
            0,
            sum(
                adjacent.GetAtomicNum() == 6
                and adjacent.GetIdx() != center
                for adjacent in neighbor.GetNeighbors()
            )
            - 1,
        )
        for neighbor in carbon_neighbors
    )
    beta_hydrogen_count = sum(
        int(neighbor.GetTotalNumHs()) for neighbor in carbon_neighbors
    )
    benzylic = any(neighbor.GetIsAromatic() for neighbor in carbon_neighbors)
    allylic = _adjacent_pi(atom, "DOUBLE")
    propargylic = _adjacent_pi(atom, "TRIPLE")
    adjacent_heteroatoms = tuple(
        sorted(
            neighbor.GetSymbol()
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() not in {1, 6}
            and neighbor.GetIdx() not in excluded
        )
    )
    ring_sizes = tuple(
        sorted(
            {
                len(ring)
                for ring in mol.GetRingInfo().AtomRings()
                if center in ring
            }
        )
    )
    contributions = tuple(
        sorted(
            (
                branch_contribution(
                    mol,
                    neighbor.GetIdx(),
                    blocked=(center, *excluded),
                    relation="alpha_carbon_branch",
                    radius=2,
                )
                for neighbor in carbon_neighbors
            ),
            key=lambda item: item.origin_atom_index,
        )
    )
    electronic = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() in excluded or neighbor.GetAtomicNum() in {1, 6}:
            continue
        electronegative = neighbor.GetSymbol() in {"F", "Cl", "Br", "I", "O", "N", "S"}
        electronic.append(
            ElectronicContribution(
                source_id=f"adjacent_heteroatom:{neighbor.GetSymbol()}",
                effect="withdrawing" if electronegative else "mixed",
                pathway="inductive",
                positional_relation="alpha",
                contribution=0.28 if electronegative else 0.0,
                atom_indices=(neighbor.GetIdx(),),
            )
        )
    context = AlkylContextDescriptor(
        context_kind="alkyl",
        carbon_substitution=carbon_substitution(atom),
        alpha_carbon_neighbor_count=len(carbon_neighbors),
        alpha_branched=len(carbon_neighbors) >= 2,
        beta_branch_count=beta_branch_count,
        beta_hydrogen_count=beta_hydrogen_count,
        cyclic=bool(ring_sizes),
        ring_sizes=ring_sizes,
        benzylic=benzylic,
        allylic=allylic,
        propargylic=propargylic,
        adjacent_heteroatoms=adjacent_heteroatoms,
    )
    return context, contributions, tuple(
        sorted(electronic, key=lambda item: (item.source_id, item.atom_indices))
    )


__all__ = ["build_alkyl_context"]
