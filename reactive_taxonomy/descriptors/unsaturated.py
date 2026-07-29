"""Alkenyl and alkynyl context descriptors."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from .common import branch_contribution
from .models import (
    AlkenylContextDescriptor,
    AlkynylContextDescriptor,
    ElectronicContribution,
    StericContribution,
)


def _pi_partner(mol: Any, center: int, order: str) -> Any | None:
    atom = mol.GetAtomWithIdx(center)
    for bond in atom.GetBonds():
        if str(bond.GetBondType()).upper() != order:
            continue
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 6:
            return other
    return None


def _endpoint_substitution(atom: Any, other_index: int) -> int:
    return sum(
        neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() != other_index
        for neighbor in atom.GetNeighbors()
    )


def _conjugation(atom_1: Any, atom_2: Any) -> str:
    endpoints = {atom_1.GetIdx(), atom_2.GetIdx()}
    if any(
        neighbor.GetIsAromatic()
        for atom in (atom_1, atom_2)
        for neighbor in atom.GetNeighbors()
        if neighbor.GetIdx() not in endpoints
    ):
        return "aryl_conjugated"
    if any(
        str(bond.GetBondType()).upper() in {"DOUBLE", "TRIPLE"}
        for atom in (atom_1, atom_2)
        for neighbor in atom.GetNeighbors()
        if neighbor.GetIdx() not in endpoints
        for bond in neighbor.GetBonds()
        if bond.GetOtherAtomIdx(neighbor.GetIdx()) not in endpoints
    ):
        return "extended_pi"
    return "isolated"


def _steric_branches(
    mol: Any, endpoints: Tuple[int, int], excluded: Iterable[int]
) -> Tuple[StericContribution, ...]:
    excluded_set = {int(value) for value in excluded} | set(endpoints)
    values = []
    for endpoint in endpoints:
        for neighbor in mol.GetAtomWithIdx(endpoint).GetNeighbors():
            if (
                neighbor.GetAtomicNum() <= 1
                or neighbor.GetIdx() in excluded_set
            ):
                continue
            values.append(
                branch_contribution(
                    mol,
                    neighbor.GetIdx(),
                    blocked=excluded_set,
                    relation="pi_endpoint_substituent",
                    radius=2,
                )
            )
    return tuple(
        sorted(values, key=lambda item: (item.origin_atom_index, item.relation))
    )


def build_alkenyl_context(
    mol: Any,
    center: int,
    *,
    excluded_atoms: Iterable[int] = (),
) -> Tuple[
    AlkenylContextDescriptor,
    Tuple[StericContribution, ...],
    Tuple[ElectronicContribution, ...],
]:
    """Build a local non-aromatic alkene context."""
    atom = mol.GetAtomWithIdx(center)
    partner = _pi_partner(mol, center, "DOUBLE")
    if partner is None:
        raise ValueError("alkenyl context requires a carbon-carbon double bond")
    endpoints = (center, partner.GetIdx())
    substitutions = (
        _endpoint_substitution(atom, partner.GetIdx()),
        _endpoint_substitution(partner, center),
    )
    total = sum(substitutions)
    alkene_class = {
        0: "unsubstituted",
        1: "monosubstituted",
        2: "disubstituted",
        3: "trisubstituted",
    }.get(total, "tetrasubstituted")
    bond = mol.GetBondBetweenAtoms(*endpoints)
    stereo = str(bond.GetStereo()) if bond is not None else None
    if stereo in {"STEREONONE", "STEREOANY"}:
        stereo = None
    ring_sizes = tuple(
        sorted(
            {
                len(ring)
                for ring in mol.GetRingInfo().AtomRings()
                if set(endpoints) <= set(ring)
            }
        )
    )
    context = AlkenylContextDescriptor(
        context_kind="alkenyl",
        endpoint_substitution=substitutions,
        alkene_class=alkene_class,
        stereochemistry=stereo,
        cyclic=bool(ring_sizes),
        ring_size=min(ring_sizes) if ring_sizes else None,
        conjugation_class=_conjugation(atom, partner),
        allylic_branch_count=sum(
            max(0, contribution.branch_count)
            for contribution in _steric_branches(
                mol, endpoints, excluded_atoms
            )
        ),
    )
    return context, _steric_branches(mol, endpoints, excluded_atoms), ()


def build_alkynyl_context(
    mol: Any,
    center: int,
    *,
    excluded_atoms: Iterable[int] = (),
) -> Tuple[
    AlkynylContextDescriptor,
    Tuple[StericContribution, ...],
    Tuple[ElectronicContribution, ...],
]:
    """Build a local non-aromatic alkyne context."""
    atom = mol.GetAtomWithIdx(center)
    partner = _pi_partner(mol, center, "TRIPLE")
    if partner is None:
        raise ValueError("alkynyl context requires a carbon-carbon triple bond")
    endpoints = (center, partner.GetIdx())
    substitutions = (
        _endpoint_substitution(atom, partner.GetIdx()),
        _endpoint_substitution(partner, center),
    )
    steric = _steric_branches(mol, endpoints, excluded_atoms)
    context = AlkynylContextDescriptor(
        context_kind="alkynyl",
        terminal=bool(atom.GetTotalNumHs() or partner.GetTotalNumHs()),
        endpoint_substitution=substitutions,
        conjugation_class=_conjugation(atom, partner),
        propargylic_branch_count=sum(item.branch_count for item in steric),
    )
    return context, steric, ()


__all__ = ["build_alkenyl_context", "build_alkynyl_context"]
