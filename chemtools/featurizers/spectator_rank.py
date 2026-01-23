"""
Spectator group ranking utilities.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Iterable, List, Set

_HETEROCYCLE_WEIGHT = 5.0
_AMINE_WEIGHT = 4.0
_ACID_WEIGHT = 3.0
_ALCOHOL_WEIGHT = 2.0
_THIOL_WEIGHT = 1.5
_DEFAULT_WEIGHT = 1.0

_AMINE_GROUPS = {
    "NH2",
    "NHR",
    "NR2",
    "NHAr",
    "NAr2",
    "AromN",
    "AromN-H",
    "Hydrazine",
}

_ACID_GROUPS = {
    "CO2H",
    "SO3H",
    "PO3H2",
    "PO3H",
    "PO2H2",
}

_ALCOHOL_GROUPS = {"OH"}
_THIOL_GROUPS = {"SH"}


@lru_cache(maxsize=1)
def _load_scaffold_motif_ids() -> Set[str]:
    path = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    motifs = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = str(entry.get("id") or "").strip()
        if motif_id:
            motifs.add(motif_id)
    return motifs


def spectator_group_weight(group_id: str) -> float:
    text = str(group_id or "").strip()
    if not text:
        return 0.0
    if text in _load_scaffold_motif_ids():
        return _HETEROCYCLE_WEIGHT
    if text in _AMINE_GROUPS:
        return _AMINE_WEIGHT
    if text in _ACID_GROUPS:
        return _ACID_WEIGHT
    if text in _ALCOHOL_GROUPS:
        return _ALCOHOL_WEIGHT
    if text in _THIOL_GROUPS:
        return _THIOL_WEIGHT
    return _DEFAULT_WEIGHT


def rank_spectator_groups(groups: Iterable[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for item in groups:
        text = str(item or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    ordered.sort(key=lambda g: (-spectator_group_weight(g), g))
    return ordered


def weighted_spectator_similarity(query_groups: Set[str], row_groups: Set[str]) -> float:
    if not query_groups and not row_groups:
        return 1.0
    if not query_groups:
        return 0.7
    if not row_groups:
        return 0.3
    intersection = query_groups & row_groups
    union = query_groups | row_groups
    if not union:
        return 0.5
    weight_union = sum(spectator_group_weight(group) for group in union)
    if weight_union <= 0:
        return 0.0
    weight_intersection = sum(spectator_group_weight(group) for group in intersection)
    return weight_intersection / weight_union


__all__ = [
    "rank_spectator_groups",
    "spectator_group_weight",
    "weighted_spectator_similarity",
]
