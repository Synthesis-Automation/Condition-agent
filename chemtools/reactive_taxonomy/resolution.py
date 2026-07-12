"""Generic typed candidate deduplication and taxonomy-directed ownership."""

from __future__ import annotations

from dataclasses import replace
from typing import Dict, Iterable, List, Set, Tuple

from .models import SiteCandidate
from .patterns import load_handle_patterns


def _pattern_index() -> Dict[str, dict]:
    return {str(item["id"]): item for item in load_handle_patterns()}


def _candidate_key(candidate: SiteCandidate) -> Tuple[object, ...]:
    return (
        candidate.site_type, candidate.topology, candidate.atom_indices,
        candidate.bond_indices, candidate.canonical_signature,
    )


def _priority(candidate: SiteCandidate, definitions: Dict[str, dict]) -> int:
    return max((int(definitions.get(pattern_id, {}).get("priority", 0)) for pattern_id in candidate.matched_patterns), default=0)


def resolve_candidates(raw_sites: Iterable[SiteCandidate]) -> List[SiteCandidate]:
    """Resolve exact duplicates and taxonomy-declared ownership conflicts."""
    definitions = _pattern_index()
    deduped: Dict[Tuple[object, ...], SiteCandidate] = {}
    for candidate in raw_sites:
        pattern_ids = tuple(sorted(
            set(candidate.matched_patterns),
            key=lambda pattern_id: (-int(definitions.get(pattern_id, {}).get("priority", 0)), pattern_id),
        ))
        details = dict(candidate.details)
        if pattern_ids:
            details["matched_pattern"] = pattern_ids[0]
            details["alternative_patterns"] = list(pattern_ids[1:])
        normalized = replace(candidate, matched_patterns=pattern_ids, details=details)
        key = _candidate_key(normalized)
        existing = deduped.get(key)
        if existing is None or _priority(normalized, definitions) > _priority(existing, definitions):
            deduped[key] = normalized

    candidates = list(deduped.values())
    suppressed: Set[int] = set()
    for owner_index, owner in enumerate(candidates):
        for pattern_id in owner.matched_patterns:
            for rule in definitions.get(pattern_id, {}).get("suppresses") or []:
                owned = set(owner.atom_roles.get(str(rule.get("owned_role") or ""), ()))
                if not owned:
                    continue
                target_type = str(rule.get("site_type") or "")
                for target_index, target in enumerate(candidates):
                    if target_index == owner_index or target.site_type != target_type:
                        continue
                    target_center = set(target.atom_roles.get("center", target.atom_indices))
                    if target_center & owned:
                        suppressed.add(target_index)

    resolved = [candidate for index, candidate in enumerate(candidates) if index not in suppressed]
    return sorted(resolved, key=lambda candidate: (candidate.atom_indices, candidate.site_type, candidate.canonical_signature))


__all__ = ["resolve_candidates"]
