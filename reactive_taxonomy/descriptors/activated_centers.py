"""Descriptors for acyl, sulfonyl, and phosphoryl reactive centers."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from .common import branch_contribution
from .models import (
    ActivatedCenterContextDescriptor,
    ElectronicContribution,
    StericContribution,
)


def _double_oxo_count(atom: Any) -> int:
    return sum(
        neighbor.GetAtomicNum() == 8
        and str(atom.GetOwningMol().GetBondBetweenAtoms(
            atom.GetIdx(), neighbor.GetIdx()
        ).GetBondType()).upper()
        == "DOUBLE"
        for neighbor in atom.GetNeighbors()
    )


def activated_context_kind(atom: Any) -> str | None:
    """Return an activated-center context kind from local atom properties."""
    element = atom.GetSymbol()
    oxo = _double_oxo_count(atom)
    if element == "C" and oxo >= 1:
        return "acyl"
    if element == "S" and oxo >= 2:
        return "sulfonyl"
    if element == "P" and oxo >= 1:
        return "phosphoryl"
    return None


def _acyl_class(atom: Any) -> str:
    neighbors = tuple(
        neighbor
        for neighbor in atom.GetNeighbors()
        if not (
            neighbor.GetAtomicNum() == 8
            and str(atom.GetOwningMol().GetBondBetweenAtoms(
                atom.GetIdx(), neighbor.GetIdx()
            ).GetBondType()).upper()
            == "DOUBLE"
        )
    )
    symbols = {neighbor.GetSymbol() for neighbor in neighbors}
    if symbols & {"F", "Cl", "Br", "I"}:
        return "acyl_halide"
    if any(
        neighbor.GetSymbol() == "O"
        and any(
            adjacent.GetAtomicNum() == 6
            and any(
                str(bond.GetBondType()).upper() == "DOUBLE"
                and bond.GetOtherAtom(adjacent).GetAtomicNum() == 8
                for bond in adjacent.GetBonds()
            )
            for adjacent in neighbor.GetNeighbors()
            if adjacent.GetIdx() != atom.GetIdx()
        )
        for neighbor in neighbors
    ):
        return "anhydride"
    if "N" in symbols:
        return "amide"
    if "O" in symbols:
        return "ester_or_acid"
    carbon_count = sum(neighbor.GetAtomicNum() == 6 for neighbor in neighbors)
    h_count = int(atom.GetTotalNumHs())
    if h_count:
        return "aldehyde"
    if carbon_count >= 2:
        return "ketone"
    return "carbonyl"


def build_activated_center_context(
    mol: Any,
    center: int,
    *,
    excluded_atoms: Iterable[int] = (),
) -> Tuple[
    ActivatedCenterContextDescriptor,
    Tuple[StericContribution, ...],
    Tuple[ElectronicContribution, ...],
]:
    """Build an activated-center profile with approach and resonance evidence."""
    atom = mol.GetAtomWithIdx(center)
    kind = activated_context_kind(atom)
    if kind is None:
        raise ValueError("center is not acyl, sulfonyl, or phosphoryl")
    excluded = {int(value) for value in excluded_atoms}
    retained = tuple(
        neighbor
        for neighbor in atom.GetNeighbors()
        if neighbor.GetIdx() not in excluded
        and not (
            neighbor.GetAtomicNum() == 8
            and str(mol.GetBondBetweenAtoms(
                center, neighbor.GetIdx()
            ).GetBondType()).upper()
            == "DOUBLE"
        )
    )
    group_classes = tuple(
        sorted(
            "aryl"
            if neighbor.GetIsAromatic()
            else "alkyl"
            if neighbor.GetAtomicNum() == 6
            else neighbor.GetSymbol()
            for neighbor in retained
        )
    )
    heteroatom_substitution = tuple(
        sorted(
            neighbor.GetSymbol()
            for neighbor in retained
            if neighbor.GetAtomicNum() not in {1, 6}
        )
    )
    steric = tuple(
        sorted(
            (
                branch_contribution(
                    mol,
                    neighbor.GetIdx(),
                    blocked=(center, *excluded),
                    relation="center_substituent",
                    radius=2,
                )
                for neighbor in retained
                if neighbor.GetAtomicNum() > 1
            ),
            key=lambda item: item.origin_atom_index,
        )
    )
    if kind == "acyl":
        center_class = _acyl_class(atom)
        activation = {
            "acyl_halide": 0.75,
            "anhydride": 0.6,
            "aldehyde": 0.4,
            "ketone": 0.25,
            "ester_or_acid": 0.05,
            "amide": -0.35,
        }.get(center_class, 0.0)
    else:
        center_class = kind
        activation = 0.5 if kind == "sulfonyl" else 0.35
    electronic = (
        ElectronicContribution(
            source_id=f"activated_center:{center_class}",
            effect="withdrawing" if activation >= 0 else "donating",
            pathway="resonance",
            positional_relation="center",
            contribution=activation,
            atom_indices=(center,),
        ),
    )
    alpha_carbons = [
        neighbor
        for neighbor in retained
        if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic()
    ]
    enolizable = (
        any(neighbor.GetTotalNumHs() > 0 for neighbor in alpha_carbons)
        if kind == "acyl"
        else None
    )
    context = ActivatedCenterContextDescriptor(
        context_kind=kind,  # type: ignore[arg-type]
        center_class=center_class,
        attached_group_classes=group_classes,
        conjugation_class=(
            "aryl_conjugated"
            if any(neighbor.GetIsAromatic() for neighbor in retained)
            else "not_aryl_conjugated"
        ),
        heteroatom_substitution=heteroatom_substitution,
        enolizable=enolizable,
    )
    return context, steric, electronic


__all__ = ["activated_context_kind", "build_activated_center_context"]
