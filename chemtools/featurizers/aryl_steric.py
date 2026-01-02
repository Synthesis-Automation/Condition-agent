"""
Aryl steric analysis around the ipso carbon using ortho bulk heuristics.
"""

from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Optional, Set

_BULK_CAP_PER_SUB = 12
_BULK_CAP_TOTAL = 20
_GROUP_BULK_CAP = 20
_STERICS_METHOD = "ortho_bulk_scaffold_v1"


def analyze_aryl_steric(mol: Any, hit: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute steric bulk around the ipso carbon for Ar-* motifs.
    """
    ipso = hit.get("a_atom_idx")
    x_root = hit.get("b_atom_idx")
    if ipso is None or x_root is None:
        return {
            "scaffold_score_0_10": 0.0,
            "score_0_10": 0.0,
            "group_bulk_0_10": 0.0,
            "description": "no steric",
            "method": _STERICS_METHOD,
            "ortho": [],
            "ortho_position_count": 0,
            "ortho_substituent_count": 0,
        }

    ring_atoms = _find_aryl_ring(mol, ipso)
    if not ring_atoms:
        return {
            "scaffold_score_0_10": 0.0,
            "score_0_10": 0.0,
            "group_bulk_0_10": 0.0,
            "description": "no steric",
            "method": _STERICS_METHOD,
            "ortho": [],
            "ortho_position_count": 0,
            "ortho_substituent_count": 0,
        }

    ortho_atoms = [
        nbr.GetIdx()
        for nbr in mol.GetAtomWithIdx(ipso).GetNeighbors()
        if nbr.GetIdx() in ring_atoms
    ]
    ring_info = mol.GetRingInfo()
    ortho_entries: List[Dict[str, Any]] = []
    for o_idx in ortho_atoms:
        ortho_entries.append(_ortho_bulk(mol, ring_atoms, o_idx, ring_info, x_root))

    bulk_total = sum(entry["bulk"] for entry in ortho_entries)
    if bulk_total > _BULK_CAP_TOTAL:
        bulk_total = _BULK_CAP_TOTAL
    scaffold_score = round(10.0 * bulk_total / _BULK_CAP_TOTAL, 1)
    group_bulk = _group_bulk_score(mol, x_root, ipso)
    ortho_substituent_count = sum(1 for entry in ortho_entries if entry.get("bulk", 0) > 0)

    if scaffold_score < 1.0:
        desc = "no steric"
    elif scaffold_score < 3.0:
        desc = "moderately steric"
    else:
        desc = "highly steric"

    return {
        "scaffold_score_0_10": scaffold_score,
        "score_0_10": scaffold_score,
        "group_bulk_0_10": group_bulk,
        "description": desc,
        "method": _STERICS_METHOD,
        "ortho": ortho_entries,
        "ortho_position_count": len(ortho_entries),
        "ortho_substituent_count": ortho_substituent_count,
    }


def _find_aryl_ring(mol: Any, ipso: int) -> Optional[Set[int]]:
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    candidates = []
    for ring in rings:
        if ipso not in ring:
            continue
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        candidates.append(ring)
    if not candidates:
        return None

    six_membered = [ring for ring in candidates if len(ring) == 6]
    target = min(six_membered or candidates, key=len)
    return set(target)


def _ortho_bulk(
    mol: Any,
    ring_atoms: Set[int],
    ortho_idx: int,
    ring_info: Any,
    x_root: int,
) -> Dict[str, Any]:
    ortho_atom = mol.GetAtomWithIdx(ortho_idx)
    subs = [
        nbr
        for nbr in ortho_atom.GetNeighbors()
        if nbr.GetIdx() not in ring_atoms and nbr.GetIdx() != x_root
    ]
    if not subs:
        return {
            "ring_atom": ortho_idx,
            "bulk": 0,
            "heavy_atoms": 0,
            "has_ring": False,
            "branching": 0,
        }

    total_bulk = 0
    total_heavy = 0
    has_ring = False
    branching = 0

    for sub in subs:
        frag = _collect_fragment(mol, sub.GetIdx(), ring_atoms)
        heavy_atoms = sum(1 for idx in frag if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
        frag_has_ring = any(mol.GetAtomWithIdx(idx).IsInRing() for idx in frag)
        frag_branching = 1 if _is_branching(mol, sub.GetIdx(), ortho_idx) else 0

        bulk = heavy_atoms
        if frag_has_ring:
            bulk += 2
        if frag_branching:
            bulk += 1
        if bulk > _BULK_CAP_PER_SUB:
            bulk = _BULK_CAP_PER_SUB

        total_bulk += bulk
        total_heavy += heavy_atoms
        has_ring = has_ring or frag_has_ring
        branching = max(branching, frag_branching)

    return {
        "ring_atom": ortho_idx,
        "bulk": total_bulk,
        "heavy_atoms": total_heavy,
        "has_ring": has_ring,
        "branching": branching,
    }


def _collect_fragment(mol: Any, start_idx: int, ring_atoms: Set[int]) -> Set[int]:
    visited: Set[int] = {start_idx}
    queue: deque[int] = deque([start_idx])
    while queue:
        idx = queue.popleft()
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx in ring_atoms or n_idx in visited:
                continue
            visited.add(n_idx)
            queue.append(n_idx)
    return visited


def _is_branching(mol: Any, start_idx: int, ring_neighbor_idx: int) -> bool:
    atom = mol.GetAtomWithIdx(start_idx)
    heavy_neighbors = [
        nbr
        for nbr in atom.GetNeighbors()
        if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != ring_neighbor_idx
    ]
    return len(heavy_neighbors) > 2


def _group_bulk_score(mol: Any, x_root: int, ipso: int) -> float:
    if x_root == ipso:
        return 0.0
    fragment = _collect_group_fragment(mol, x_root, ipso)
    heavy_atoms = sum(1 for idx in fragment if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
    has_ring = any(mol.GetAtomWithIdx(idx).IsInRing() for idx in fragment)
    root_atom = mol.GetAtomWithIdx(x_root)
    root_heavy_neighbors = [
        nbr for nbr in root_atom.GetNeighbors()
        if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != ipso
    ]
    branching = 1 if len(root_heavy_neighbors) >= 2 else 0

    bulk = heavy_atoms + (2 if has_ring else 0) + branching
    bulk = min(bulk, _GROUP_BULK_CAP)
    return round(10.0 * bulk / _GROUP_BULK_CAP, 1)


def _collect_group_fragment(mol: Any, start_idx: int, ipso: int) -> Set[int]:
    visited: Set[int] = {start_idx}
    queue: deque[int] = deque([start_idx])
    while queue:
        idx = queue.popleft()
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx == ipso or n_idx in visited:
                continue
            visited.add(n_idx)
            queue.append(n_idx)
    return visited
