"""Exact bond differences for reaction SMILES with supplied atom maps."""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

from .chemistry.rdkit_utils import parse_smiles


def _mapped_bonds(side: str) -> Dict[Tuple[int, int], str]:
    bonds: Dict[Tuple[int, int], str] = {}
    for token in (part for part in side.split(".") if part):
        mol = parse_smiles(token)
        if mol is None: continue
        for bond in mol.GetBonds():
            left, right = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if not left or not right: continue
            bonds[tuple(sorted((int(left), int(right))))] = str(bond.GetBondType())
    return bonds


def supplied_map_bond_changes(reaction_smiles: str) -> List[Dict[str, Any]]:
    if ">>" in reaction_smiles: left, right = reaction_smiles.split(">>", 1)
    else:
        parts = reaction_smiles.split(">")
        if len(parts) != 3: return []
        left, right = parts[0], parts[2]
    reactant, product = _mapped_bonds(left), _mapped_bonds(right)
    changes: List[Dict[str, Any]] = []
    for pair in sorted(set(reactant) | set(product)):
        old, new = reactant.get(pair), product.get(pair)
        if old == new: continue
        kind = "formed" if old is None else ("broken" if new is None else "order_changed")
        changes.append({"change_type": kind, "atom_maps": list(pair), "old_order": old, "new_order": new, "evidence": "supplied_atom_mapping"})
    return changes


__all__ = ["supplied_map_bond_changes"]
