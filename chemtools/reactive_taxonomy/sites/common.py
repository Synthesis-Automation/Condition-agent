"""Helpers shared by site detectors."""

from __future__ import annotations

from typing import Any, Iterable, List


def bond_index(mol: Any, left: int, right: int) -> int:
    bond = mol.GetBondBetweenAtoms(left, right)
    return bond.GetIdx() if bond is not None else -1


def unique_indices(values: Iterable[int]) -> List[int]:
    return sorted(set(int(value) for value in values))
