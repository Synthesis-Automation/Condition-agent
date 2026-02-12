"""
Reaction-level helpers backed by taxonomy aliases.

All reaction alias normalization is sourced from
``taxonomy/data/reaction_types.v4.0.json``.
"""

from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional, Set

from ...taxonomy import reaction_catalog as _reaction_catalog


def slugify_family(value: str) -> str:
    """Return a slug suitable for alias lookup."""
    return re.sub(r"[^0-9a-z]+", "_", value.lower()).strip("_")


def _candidate_family_labels(label: str) -> list[str]:
    stripped = label.strip()
    if not stripped:
        return []
    candidates = [
        stripped,
        slugify_family(stripped),
        stripped.replace("-", "_"),
        stripped.replace(" ", "_"),
        stripped.replace("-", " "),
    ]
    seen: set[str] = set()
    ordered: list[str] = []
    for candidate in candidates:
        text = candidate.strip()
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return ordered


def _preferred_reaction_label(reaction_id: str) -> str:
    if not reaction_id:
        return reaction_id
    preferred = _reaction_label_preferences().get(reaction_id)
    if preferred:
        return preferred
    if re.fullmatch(r"[a-z0-9_]+", reaction_id):
        return reaction_id
    reaction = _reaction_catalog.get_reaction_type(reaction_id)
    if not reaction:
        return reaction_id
    for alias in reaction.aliases:
        if isinstance(alias, str) and re.fullmatch(r"[a-z0-9_]+", alias):
            return alias
    return reaction_id


@lru_cache(maxsize=1)
def _reaction_label_preferences() -> dict[str, str]:
    try:
        from ...taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_featurizer_logic()
    if not payload:
        return {}
    mapping = payload.get("reaction_label_preferences", {}) or {}
    if not isinstance(mapping, dict):
        return {}
    out: dict[str, str] = {}
    for key, value in mapping.items():
        if not isinstance(key, str) or not isinstance(value, str):
            continue
        key_norm = key.strip()
        value_norm = value.strip()
        if key_norm and value_norm:
            out[key_norm] = value_norm
    return out


def canonical_family_label(family: Optional[str]) -> Optional[str]:
    """Resolve ``family`` to a stable canonical label."""
    if not family:
        return None
    for candidate in _candidate_family_labels(family):
        canonical = _reaction_catalog.resolve_reaction_type(candidate)
        if canonical:
            return _preferred_reaction_label(canonical)
    return None


def resolve_reaction_family(family: Optional[str]) -> Optional[str]:
    """Resolve arbitrary family label/alias to canonical taxonomy label."""
    return canonical_family_label(family)


def _canonical_label_set(labels: list[str]) -> Set[str]:
    values: Set[str] = set()
    for label in labels:
        canonical = canonical_family_label(label)
        if canonical:
            values.add(canonical)
    return values


@lru_cache(maxsize=1)
def _canonical_cn_family() -> str:
    for candidate in ("c_n_cross_coupling", "C_N_Coupling"):
        canonical = canonical_family_label(candidate)
        if canonical:
            return canonical
    return "c_n_cross_coupling"


CN_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(
    ["c_n_cross_coupling", "C_N_Coupling", "snar_cn", "Snar_cn"]
)
CO_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["c_o_coupling", "C_O_Coupling"])
CS_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["c_s_coupling", "C_S_Coupling"])
AMIDE_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(
    ["amide_formation", "Amide_formation", "amide_coupling"]
)
SONOGASHIRA_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["sonogashira", "Sonogashira"])
SUZUKI_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["suzuki_miyaura", "Suzuki_miyaura"])
NEGISHI_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["negishi", "Negishi"])
STILLE_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["stille", "Stille"])
KUMADA_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["kumada", "Kumada"])
HIYAMA_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["hiyama", "Hiyama"])
HECK_FAMILIES_CANONICAL: Set[str] = _canonical_label_set(["heck", "Heck"])
RCM_FAMILIES_CANONICAL: Set[str] = set()
ULLMANN_SPECIFIC_CANONICAL: Set[str] = set()
BUCHWALD_SPECIFIC_CANONICAL: Set[str] = set()


def apply_catalyst_override(
    family: str,
    metals: Set[str],
    *,
    is_cn_coupling: bool,
) -> str:
    """Apply catalyst-based family override for C-N coupling reactions."""
    canonical = canonical_family_label(family) or family or "Unknown"
    if not metals:
        return canonical
    if canonical in CN_FAMILIES_CANONICAL or is_cn_coupling:
        return _canonical_cn_family()
    return canonical


__all__ = [
    "CN_FAMILIES_CANONICAL",
    "CO_FAMILIES_CANONICAL",
    "CS_FAMILIES_CANONICAL",
    "AMIDE_FAMILIES_CANONICAL",
    "SONOGASHIRA_FAMILIES_CANONICAL",
    "SUZUKI_FAMILIES_CANONICAL",
    "NEGISHI_FAMILIES_CANONICAL",
    "STILLE_FAMILIES_CANONICAL",
    "KUMADA_FAMILIES_CANONICAL",
    "HIYAMA_FAMILIES_CANONICAL",
    "HECK_FAMILIES_CANONICAL",
    "RCM_FAMILIES_CANONICAL",
    "ULLMANN_SPECIFIC_CANONICAL",
    "BUCHWALD_SPECIFIC_CANONICAL",
    "slugify_family",
    "canonical_family_label",
    "resolve_reaction_family",
    "apply_catalyst_override",
]
