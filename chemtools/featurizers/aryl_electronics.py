"""
Tag-weighted aryl electronics analysis around an ipso carbon.
"""

from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Optional, Set, Tuple

from chemtools.util.smarts_cache import compile_smarts

_ELECTRONICS_METHOD = "tag_weighted_v1"
_ELECTRONIC_GROUP_IDS = ["NO2", "CN", "CHO", "CO2R", "CO2H", "CONHR"]
_WEIGHTS = {"ortho": 1.0, "meta": 0.6, "para": 1.2, "ipso": 1.0}
_STRENGTHS = {"ewg_strong": 2.0, "ewg_moderate": 1.0}


def analyze_aryl_electronics(
    mol: Any,
    hit: Dict[str, Any],
    groups_dict: Dict[str, Dict[str, Any]],
    *,
    include_ipso_group: bool = True,
    include_gasteiger: bool = False,
) -> Dict[str, Any]:
    """
    Compute tag-weighted electronics around the ipso carbon.
    """
    ipso = hit.get("a_atom_idx")
    if ipso is None:
        return {
            "score_0_10": 5.0,
            "method": _ELECTRONICS_METHOD,
            "including_group": include_ipso_group,
            "contributions": [],
        }

    ring_atoms = _find_aryl_ring(mol, ipso)
    if not ring_atoms:
        return {
            "score_0_10": 5.0,
            "method": _ELECTRONICS_METHOD,
            "including_group": include_ipso_group,
            "contributions": [],
        }

    ring_dist = _ring_distances(mol, ring_atoms, ipso)
    group_by_atom = _map_group_matches(mol, groups_dict)
    contributions: List[Dict[str, Any]] = []
    e_sum = 0.0

    for ring_atom in sorted(ring_atoms, key=lambda idx: (ring_dist.get(idx, 99), idx)):
        if ring_atom != ipso and ring_atom not in ring_dist:
            continue
        for nbr in mol.GetAtomWithIdx(ring_atom).GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in ring_atoms:
                continue
            group_id = group_by_atom.get(nbr_idx)
            if not group_id:
                continue
            if ring_atom == ipso and not include_ipso_group:
                continue

            pos = "ipso" if ring_atom == ipso else _position_for_distance(ring_dist[ring_atom])
            if pos is None:
                continue
            strength = _strength_from_tags(groups_dict.get(group_id, {}).get("tags") or [])
            weight = _WEIGHTS.get(pos, 0.0)
            if strength == 0.0 or weight == 0.0:
                continue
            term = round(weight * strength, 3)
            e_sum += term
            contributions.append(
                {
                    "pos": pos,
                    "group_id": group_id,
                    "strength": strength,
                    "weight": weight,
                    "term": term,
                }
            )

    score = _clamp(round(5.0 + e_sum, 1), 0.0, 10.0)
    result: Dict[str, Any] = {
        "score_0_10": score,
        "method": _ELECTRONICS_METHOD,
        "including_group": include_ipso_group,
        "contributions": contributions,
    }

    if include_gasteiger:
        gasteiger = _gasteiger_charge(mol, ipso)
        if gasteiger is not None:
            result["optional"] = {"gasteiger": {"q_ipso": gasteiger}}

    return result


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


def _position_for_distance(distance: int) -> Optional[str]:
    if distance == 1:
        return "ortho"
    if distance == 2:
        return "meta"
    if distance == 3:
        return "para"
    return None


def _map_group_matches(mol: Any, groups_dict: Dict[str, Dict[str, Any]]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for group_id in _ELECTRONIC_GROUP_IDS:
        group = groups_dict.get(group_id)
        if not group:
            continue
        smarts = group.get("smarts", "")
        if not smarts:
            continue
        pattern = compile_smarts(smarts, validate=False)
        if pattern is None:
            continue
        map_idx = _query_map_atom(pattern, map_num=2)
        if map_idx is None:
            continue
        matches = mol.GetSubstructMatches(pattern, uniquify=True)
        for match in matches:
            atom_idx = match[map_idx]
            if atom_idx not in mapping:
                mapping[atom_idx] = group_id
    return mapping


def _query_map_atom(pattern: Any, *, map_num: int) -> Optional[int]:
    for atom in pattern.GetAtoms():
        if atom.GetAtomMapNum() == map_num:
            return atom.GetIdx()
    return None


def _strength_from_tags(tags: List[str]) -> float:
    for tag in tags:
        if tag in _STRENGTHS:
            return _STRENGTHS[tag]
    return 0.0


def _gasteiger_charge(mol: Any, ipso: int) -> Optional[float]:
    try:
        from rdkit.Chem import rdPartialCharges
    except Exception:
        return None
    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
        charge = mol.GetAtomWithIdx(ipso).GetProp("_GasteigerCharge")
        return float(charge)
    except Exception:
        return None


def _clamp(value: float, min_value: float, max_value: float) -> float:
    return max(min_value, min(max_value, value))
