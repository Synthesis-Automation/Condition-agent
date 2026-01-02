"""
Simplified nearby functional group listing relative to a reaction center.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

def analyze_nearby_groups(
    mol: Any,
    hit: Dict[str, Any],
    all_motifs: List[Dict[str, Any]],
    groups_dict: Dict[str, Dict[str, Any]],
    compound_map: Optional[Dict[str, Any]] = None,
) -> List[str]:
    """
    Identify functional groups near the reaction center (simple labels).
    """
    center_ipso = hit.get("a_atom_idx")
    center_fg = hit.get("b_atom_idx")
    if center_ipso is None:
        return []

    seen: Set[str] = set()
    labels: List[str] = []

    for other in all_motifs:
        # Skip the current reaction center itself
        if other.get("a_atom_idx") == center_ipso and other.get("b_atom_idx") == center_fg:
            continue
            
        other_id = other["compound_id"]
        other_a = other.get("a_atom_idx")
        if other_a is None:
            continue
        
        group_id = _resolve_group_id(other_id, compound_map)
        group_info = groups_dict.get(group_id, {}) if group_id else {}
        label = group_info.get("name") or group_id or other_id

        if group_id == "H" or label == "-H":
            continue

        if label not in seen:
            seen.add(label)
            labels.append(label)

    return sorted(labels)


def _resolve_group_id(compound_id: str, compound_map: Optional[Dict[str, Any]]) -> Optional[str]:
    if compound_map and compound_id in compound_map:
        entry = compound_map[compound_id]
        group_b = getattr(entry, "group_b", None)
        if group_b:
            return str(group_b)
    if "-" in compound_id:
        return compound_id.split("-")[-1]
    return None
