"""
Shared ring utilities for featurizers.
"""

from __future__ import annotations

from typing import Any, Optional, Set


def find_aryl_ring(mol: Any, ipso: int) -> Optional[Set[int]]:
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
