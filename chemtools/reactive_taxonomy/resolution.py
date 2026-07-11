"""Generic candidate deduplication and taxonomy-directed site ownership."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Set, Tuple

from .patterns import load_handle_patterns


def _pattern_index() -> Dict[str, Dict[str, Any]]:
    return {str(item["id"]): item for item in load_handle_patterns()}


def _candidate_key(site_type: str, item: Dict[str, Any]) -> Tuple[Any, ...]:
    return (
        site_type,
        item.get("topology"),
        tuple(item.get("atom_indices") or []),
        tuple(item.get("bond_indices") or []),
        item.get("signature"),
    )


def _priority(item: Dict[str, Any], definitions: Dict[str, Dict[str, Any]]) -> int:
    return max(
        (int(definitions.get(pattern_id, {}).get("priority", 0)) for pattern_id in item.get("matched_patterns") or []),
        default=0,
    )


def _owned_atoms(item: Dict[str, Any], role: str) -> Set[int]:
    details = item.get("details") or {}
    role_fields = {
        "handle": "handle_atom_indices",
        "center": "center_atom_index",
        "anchor": "anchor_atom_index",
        "connector": "connector_atom_index",
    }
    value = details.get(role_fields.get(role, ""))
    if value is None:
        return set()
    if isinstance(value, list):
        return {int(index) for index in value}
    return {int(value)}


def resolve_candidates(raw_sites: Iterable[Tuple[str, Dict[str, Any]]]) -> List[Tuple[str, Dict[str, Any]]]:
    """Resolve exact duplicates and taxonomy-declared ownership conflicts."""
    definitions = _pattern_index()
    deduped: Dict[Tuple[Any, ...], Tuple[str, Dict[str, Any]]] = {}
    for site_type, raw_item in raw_sites:
        item = dict(raw_item)
        pattern_ids = sorted(
            set(str(value) for value in item.get("matched_patterns") or []),
            key=lambda pattern_id: (-int(definitions.get(pattern_id, {}).get("priority", 0)), pattern_id),
        )
        item["matched_patterns"] = pattern_ids
        if pattern_ids:
            details = dict(item.get("details") or {})
            details["matched_pattern"] = pattern_ids[0]
            details["alternative_patterns"] = pattern_ids[1:]
            item["details"] = details
        key = _candidate_key(site_type, item)
        existing = deduped.get(key)
        if existing is None or _priority(item, definitions) > _priority(existing[1], definitions):
            deduped[key] = (site_type, item)

    candidates = list(deduped.values())
    suppressed: Set[int] = set()
    for owner_index, (_, owner) in enumerate(candidates):
        for pattern_id in owner.get("matched_patterns") or []:
            definition = definitions.get(pattern_id, {})
            for rule in definition.get("suppresses") or []:
                owned = _owned_atoms(owner, str(rule.get("owned_role") or ""))
                if not owned:
                    continue
                target_type = str(rule.get("site_type") or "")
                for target_index, (site_type, target) in enumerate(candidates):
                    if target_index == owner_index or site_type != target_type:
                        continue
                    target_center = _owned_atoms(target, "center") or set(target.get("atom_indices") or [])
                    if target_center & owned:
                        suppressed.add(target_index)

    resolved = [candidate for index, candidate in enumerate(candidates) if index not in suppressed]
    resolved.sort(key=lambda pair: (pair[1].get("atom_indices") or [], pair[0], pair[1].get("signature") or ""))
    return resolved


__all__ = ["resolve_candidates"]
