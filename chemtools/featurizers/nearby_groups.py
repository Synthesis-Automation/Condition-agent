"""
Analysis of nearby functional groups relative to a reaction center.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set
from rdkit import Chem
from rdkit.Chem import rdmolops

def analyze_nearby_groups(
    mol: Any,
    hit: Dict[str, Any],
    all_motifs: List[Dict[str, Any]],
    groups_dict: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Identify and characterize functional groups near the reaction center.
    """
    center_ipso = hit.get("a_atom_idx")
    center_fg = hit.get("b_atom_idx")
    if center_ipso is None:
        return []

    nearby = []
    
    # Get the ring atoms if it's an aryl center
    ring_atoms = _find_aryl_ring(mol, center_ipso)
    ring_dist = {}
    if ring_atoms and center_ipso in ring_atoms:
        ring_dist = _ring_distances(mol, ring_atoms, center_ipso)

    for other in all_motifs:
        # Skip the current reaction center itself
        if other.get("a_atom_idx") == center_ipso and other.get("b_atom_idx") == center_fg:
            continue
            
        other_id = other["compound_id"]
        other_a = other.get("a_atom_idx")
        if other_a is None:
            continue
        
        # Calculate topological distance safely
        if center_ipso == other_a:
            dist = 0
        else:
            try:
                path = rdmolops.GetShortestPath(mol, int(center_ipso), int(other_a))
                dist = len(path) - 1 if path else 99
            except Exception:
                # Fallback for complex fused systems or RDKit invariant violations
                dist = 99
        
        pos = None
        if ring_atoms and other_a in ring_atoms:
            d = ring_dist.get(other_a, 0)
            if d == 1: pos = "ortho"
            elif d == 2: pos = "meta"
            elif d == 3: pos = "para"
        elif dist == 1:
            pos = "alpha"
        elif dist == 2:
            pos = "beta"
        elif dist == 3:
            pos = "gamma"

        # Get tags from the group definition if available
        # Note: compound_id might be "Ar-Br", we might need to look up the group "Br"
        # or just use the compound tags.
        tags = []
        # Try to find the 'B' side group tags (the functional group part)
        # This assumes the compound_id follows the "A-B" pattern from the registry
        if "-" in other_id:
            fg_part = other_id.split("-")[-1]
            group_info = groups_dict.get(fg_part, {})
            tags = group_info.get("tags", [])

        nearby.append({
            "motif_id": other_id,
            "distance": dist,
            "position": pos,
            "tags": tags
        })

    return sorted(nearby, key=lambda x: x["distance"])


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


def _ring_distances(mol: Any, ring_atoms: Set[int], ipso: int) -> Dict[int, int]:
    from collections import deque
    dist: Dict[int, int] = {ipso: 0}
    queue: deque[int] = deque([ipso])
    while queue:
        idx = queue.popleft()
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            n_idx = nbr.GetIdx()
            if n_idx not in ring_atoms or n_idx in dist:
                continue
            dist[n_idx] = dist[idx] + 1
            queue.append(n_idx)
    return dist
