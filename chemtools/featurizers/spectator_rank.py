"""
Spectator group ranking utilities.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable, List, Set

_DEFAULT_WEIGHTS = {
    "default": 1.0,
    "scaffold": 5.0,
    "amine": 4.0,
    "acid": 3.0,
    "alcohol": 2.0,
    "thiol": 1.5,
}

_DEFAULT_AMINE_GROUPS = {
    "NH2",
    "NHR",
    "NR2",
    "NHAr",
    "NAr2",
    "AromN",
    "AromN-H",
    "Hydrazine",
}

_DEFAULT_ACID_GROUPS = {
    "CO2H",
    "SO3H",
    "PO3H2",
    "PO3H",
    "PO2H2",
}

_DEFAULT_ALCOHOL_GROUPS = {"OH"}
_DEFAULT_THIOL_GROUPS = {"SH"}


@lru_cache(maxsize=1)
def _load_spectator_logic() -> dict:
    try:
        from ..taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_featurizer_logic()
    if not payload:
        return {}
    section = payload.get("spectator_rank", {}) or {}
    return section if isinstance(section, dict) else {}


def _configured_weight(name: str) -> float:
    logic = _load_spectator_logic()
    weights = logic.get("weights", {}) or {}
    if not isinstance(weights, dict):
        return float(_DEFAULT_WEIGHTS[name])
    value = weights.get(name, _DEFAULT_WEIGHTS[name])
    try:
        return float(value)
    except Exception:
        return float(_DEFAULT_WEIGHTS[name])


def _configured_group_set(name: str, default: Set[str]) -> Set[str]:
    logic = _load_spectator_logic()
    group_sets = logic.get("group_sets", {}) or {}
    if not isinstance(group_sets, dict):
        return set(default)
    values = group_sets.get(name, [])
    if not isinstance(values, list):
        return set(default)
    configured = {str(v) for v in values if isinstance(v, str) and v.strip()}
    return configured or set(default)


@lru_cache(maxsize=1)
def _load_scaffold_motif_ids() -> Set[str]:
    try:
        from ..taxonomy import loader as taxonomy_loader
    except Exception:
        return set()
    payload = taxonomy_loader.load_scaffold_motifs()
    if not payload:
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
    amine_groups = _configured_group_set("amine", _DEFAULT_AMINE_GROUPS)
    acid_groups = _configured_group_set("acid", _DEFAULT_ACID_GROUPS)
    alcohol_groups = _configured_group_set("alcohol", _DEFAULT_ALCOHOL_GROUPS)
    thiol_groups = _configured_group_set("thiol", _DEFAULT_THIOL_GROUPS)
    if text in _load_scaffold_motif_ids():
        return _configured_weight("scaffold")
    if text in amine_groups:
        return _configured_weight("amine")
    if text in acid_groups:
        return _configured_weight("acid")
    if text in alcohol_groups:
        return _configured_weight("alcohol")
    if text in thiol_groups:
        return _configured_weight("thiol")
    return _configured_weight("default")


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
