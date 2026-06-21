"""
Alkyl steric analysis for R-/Bn-/Allyl-context motifs (agnostic to X).
"""

from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Optional, Set

_METHOD = "alkyl_alpha_scaffold_v1"
_SUB_BULK_CAP = 12
_RAW_CAP = 20
_GROUP_BULK_CAP = 20
_BASELINE = {"methyl": 0, "primary": 2, "secondary": 4, "tertiary": 7}


def analyze_alkyl_steric(
    mol: Any,
    hit: Dict[str, Any],
    *,
    include_details: bool = False,
) -> Dict[str, Any]:
    """
    Compute steric scaffold score around an alkyl alpha carbon.
    """
    def empty_result() -> Dict[str, Any]:
        result: Dict[str, Any] = {
            "score_0_10": 0.0,
            "group_bulk_0_10": 0.0,
            "description": "no steric",
            "method": _METHOD,
            "classification": "unknown",
            "beta_hydrogens": 0,
        }
        if include_details:
            result["details"] = {"alpha": {}, "substituents": []}
        return result

    alpha = hit.get("a_atom_idx")
    x_root = hit.get("b_atom_idx")
    if alpha is None or x_root is None:
        return empty_result()

    alpha_atom = mol.GetAtomWithIdx(alpha)
    neighbors = list(alpha_atom.GetNeighbors())
    subs = [nbr for nbr in neighbors if nbr.GetIdx() != x_root]

    carbon_subs = _count_carbon_subs(subs)
    classification = _classify_alpha(carbon_subs)
    benzylic = any(nbr.GetIsAromatic() for nbr in subs)
    allylic = _is_allylic(alpha_atom)
    has_beta_branching = _has_beta_branching(subs, alpha)
    beta_hydrogens = _count_beta_hydrogens(mol, alpha, x_root)

    substituents = []
    bulk_total = 0
    for sub in subs:
        entry = _substituent_bulk(mol, sub.GetIdx(), alpha, x_root)
        substituents.append(entry)
        bulk_total += entry["bulk"]

    raw = _BASELINE.get(classification, 0) + bulk_total
    raw = min(raw, _RAW_CAP)
    if _is_neopentyl(alpha_atom, subs):
        raw = min(raw + 3, _RAW_CAP)

    scaffold_score = round(10.0 * raw / _RAW_CAP, 1)
    group_bulk = _group_bulk_score(mol, x_root, alpha)

    if scaffold_score < 1.0:
        desc = "no steric"
    elif scaffold_score < 3.0:
        desc = "moderately steric"
    else:
        desc = "highly steric"

    result: Dict[str, Any] = {
        "score_0_10": scaffold_score,
        "group_bulk_0_10": group_bulk,
        "description": desc,
        "method": _METHOD,
        "classification": classification,
        "beta_hydrogens": beta_hydrogens,
    }
    if include_details:
        result["details"] = {
            "scaffold_score_0_10": scaffold_score,
            "alpha": {
                "degree": alpha_atom.GetDegree(),
                "carbon_subs": carbon_subs,
                "benzylic": benzylic,
                "allylic": allylic,
                "has_beta_branching": has_beta_branching,
                "beta_hydrogens": beta_hydrogens,
            },
            "substituents": substituents,
        }
    return result


def _count_carbon_subs(subs: List[Any]) -> int:
    return sum(1 for nbr in subs if nbr.GetAtomicNum() == 6)


def _classify_alpha(carbon_subs: int) -> str:
    if carbon_subs <= 0:
        return "methyl"
    if carbon_subs == 1:
        return "primary"
    if carbon_subs == 2:
        return "secondary"
    return "tertiary"


def _is_allylic(alpha_atom: Any) -> bool:
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetAtomicNum() != 6:
            continue
        for bond in nbr.GetBonds():
            if bond.GetBondTypeAsDouble() != 2.0:
                continue
            other = bond.GetOtherAtom(nbr)
            if other.GetIdx() == alpha_atom.GetIdx():
                continue
            if other.GetAtomicNum() == 6:
                return True
    return False


def _has_beta_branching(subs: List[Any], alpha_idx: int) -> bool:
    for sub in subs:
        heavy_neighbors = [
            nbr
            for nbr in sub.GetNeighbors()
            if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != alpha_idx
        ]
        if len(heavy_neighbors) >= 2:
            return True
    return False


def _count_beta_hydrogens(mol: Any, alpha_idx: int, x_root: int) -> int:
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    total = 0
    for nbr in alpha_atom.GetNeighbors():
        n_idx = nbr.GetIdx()
        if n_idx == x_root:
            continue
        if nbr.GetAtomicNum() != 6:
            continue
        if nbr.GetIsAromatic():
            continue
        explicit_h = sum(1 for h in nbr.GetNeighbors() if h.GetAtomicNum() == 1)
        total += explicit_h + nbr.GetNumImplicitHs()
    return total


def _substituent_bulk(mol: Any, sub_idx: int, alpha_idx: int, x_root: int) -> Dict[str, Any]:
    local = _local_radius_two(mol, sub_idx, alpha_idx, x_root)
    heavy_atoms_r2 = sum(1 for idx in local if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
    has_ring = any(mol.GetAtomWithIdx(idx).IsInRing() for idx in local)

    sub_atom = mol.GetAtomWithIdx(sub_idx)
    heavy_neighbors = [
        nbr for nbr in sub_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != alpha_idx
    ]
    branching = 1 if len(heavy_neighbors) >= 2 else 0

    bulk = heavy_atoms_r2 + (2 if has_ring else 0) + branching
    bulk = min(bulk, _SUB_BULK_CAP)

    return {
        "sub_atom": sub_idx,
        "bulk": bulk,
        "heavy_atoms_r2": heavy_atoms_r2,
        "branching": branching,
        "has_ring": has_ring,
    }


def _local_radius_two(mol: Any, sub_idx: int, alpha_idx: int, x_root: int) -> Set[int]:
    local: Set[int] = {sub_idx}
    sub_atom = mol.GetAtomWithIdx(sub_idx)
    for nbr in sub_atom.GetNeighbors():
        n_idx = nbr.GetIdx()
        if n_idx in {alpha_idx, x_root}:
            continue
        local.add(n_idx)
    return local


def _is_neopentyl(alpha_atom: Any, subs: List[Any]) -> bool:
    carbon_subs = [nbr for nbr in subs if nbr.GetAtomicNum() == 6]
    if len(carbon_subs) != 1:
        return False
    beta = carbon_subs[0]
    carbon_neighbors = [
        nbr
        for nbr in beta.GetNeighbors()
        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != alpha_atom.GetIdx()
    ]
    return len(carbon_neighbors) >= 3


def _group_bulk_score(mol: Any, x_root: int, alpha_idx: int) -> float:
    fragment = _collect_fragment(mol, x_root, alpha_idx)
    heavy_atoms = sum(1 for idx in fragment if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
    has_ring = any(mol.GetAtomWithIdx(idx).IsInRing() for idx in fragment)
    root_atom = mol.GetAtomWithIdx(x_root)
    root_heavy_neighbors = [
        nbr
        for nbr in root_atom.GetNeighbors()
        if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != alpha_idx
    ]
    branching = 1 if len(root_heavy_neighbors) >= 2 else 0

    bulk = heavy_atoms + (2 if has_ring else 0) + branching
    bulk = min(bulk, _GROUP_BULK_CAP)
    return round(10.0 * bulk / _GROUP_BULK_CAP, 1)


def _collect_fragment(mol: Any, start_idx: int, alpha_idx: int) -> Set[int]:
    visited: Set[int] = {start_idx}
    queue: deque[int] = deque([start_idx])
    while queue:
        idx = queue.popleft()
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx == alpha_idx or n_idx in visited:
                continue
            visited.add(n_idx)
            queue.append(n_idx)
    return visited
