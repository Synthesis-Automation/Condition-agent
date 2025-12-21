"""
Reaction-level helpers that interface with the taxonomy registry.

This module centralises the family alias logic that was previously embedded in
``chemtools.router`` so that other components can consume the canonical
reaction identifiers without depending on the routing heuristics.
"""

from __future__ import annotations

import re
from typing import Optional, Dict, Set

from ..taxonomy import reaction_catalog as _reaction_catalog

FAMILY_ALIAS_OVERRIDES: Dict[str, str] = {
    # C-N Coupling reactions (legacy -> v2)
    "C_N_Coupling": "c_n_cross_coupling",
    "Buchwald_CN": "c_n_cross_coupling",
    "Buchwald-Hartwig": "c_n_cross_coupling",
    "Ullmann_CN": "c_n_cross_coupling",
    "Chan_Lam_CN": "c_n_cross_coupling",
    "chan_lam_cn": "c_n_cross_coupling",
    "cn_coupling": "c_n_cross_coupling",
    "buchwald_hartwig_c_n": "c_n_cross_coupling",
    "ullmann_cn": "c_n_cross_coupling",
    
    # C-C Coupling reactions
    "Suzuki_CC": "suzuki_miyaura",
    "Suzuki": "suzuki_miyaura",
    "Suzuki-Miyaura": "suzuki_miyaura",
    "Sonogashira_CC": "sonogashira",
    "Sonogashira": "sonogashira",
    "Heck": "heck",
    
    # Amide Coupling
    "Amide_Coupling": "amide_coupling",
    "Amide_formation": "amide_coupling",
    "amide_formation": "amide_coupling",
    
    # SNAr reactions (legacy -> v2)
    "S_NAr": "snar_cn",
    "SNAr": "snar_cn",
    "snar": "snar_cn",
    "s_nar": "snar_cn",
    "Aromatic_Nucleophilic_Substitution": "snar_cn",
    "aromatic_nucleophilic_substitution": "snar_cn",
    
    # Esterification
    "Esterification": "esterification",
    "Ester_formation": "esterification",
}

CN_FAMILIES_CANONICAL: Set[str] = {"c_n_cross_coupling", "snar_cn"}
CO_FAMILIES_CANONICAL: Set[str] = set()
CS_FAMILIES_CANONICAL: Set[str] = set()
AMIDE_FAMILIES_CANONICAL: Set[str] = {"amide_coupling"}
SONOGASHIRA_FAMILIES_CANONICAL: Set[str] = {"sonogashira"}
SUZUKI_FAMILIES_CANONICAL: Set[str] = {"suzuki_miyaura"}
NEGISHI_FAMILIES_CANONICAL: Set[str] = set()
STILLE_FAMILIES_CANONICAL: Set[str] = set()
HECK_FAMILIES_CANONICAL: Set[str] = {"heck"}
RCM_FAMILIES_CANONICAL: Set[str] = set()
ULLMANN_SPECIFIC_CANONICAL: Set[str] = set()
BUCHWALD_SPECIFIC_CANONICAL: Set[str] = set()


def slugify_family(value: str) -> str:
    """Return a slug suitable for alias lookup."""
    return re.sub(r"[^0-9a-z]+", "_", value.lower()).strip("_")


def canonical_family_label(family: Optional[str]) -> Optional[str]:
    """Resolve ``family`` to a canonical taxonomy identifier when possible."""
    if not family:
        return None
    family = family.strip()
    if not family:
        return None

    canonical = _reaction_catalog.resolve_reaction_type(family)
    if canonical:
        return canonical
    override = _resolve_family_override(family)
    if override:
        return override
    slug = slugify_family(family)
    if slug and slug != family:
        canonical = _reaction_catalog.resolve_reaction_type(slug)
        if canonical:
            return canonical
        override = _resolve_family_override(slug)
        if override:
            return override
    return None


def resolve_reaction_family(family: Optional[str]) -> Optional[str]:
    """Resolve arbitrary family label/alias to canonical taxonomy ID (or None)."""
    return canonical_family_label(family)


def apply_catalyst_override(
    family: str,
    metals: Set[str],
    *,
    is_cn_coupling: bool,
) -> str:
    """
    Apply catalyst-based family override for C-N coupling reactions using canonical IDs.
    """
    canonical = family or "Unknown"
    if not metals:
        return canonical
    if canonical in CN_FAMILIES_CANONICAL or is_cn_coupling:
        return "c_n_cross_coupling"
    return canonical


def _resolve_family_override(family: str) -> Optional[str]:
    for candidate in (family, family.lower(), slugify_family(family)):
        if not candidate:
            continue
        override = FAMILY_ALIAS_OVERRIDES.get(candidate)
        if override and _reaction_catalog.get_reaction_type(override):
            return override
    return None


__all__ = [
    "FAMILY_ALIAS_OVERRIDES",
    "CN_FAMILIES_CANONICAL",
    "CO_FAMILIES_CANONICAL",
    "CS_FAMILIES_CANONICAL",
    "AMIDE_FAMILIES_CANONICAL",
    "SONOGASHIRA_FAMILIES_CANONICAL",
    "SUZUKI_FAMILIES_CANONICAL",
    "NEGISHI_FAMILIES_CANONICAL",
    "STILLE_FAMILIES_CANONICAL",
    "HECK_FAMILIES_CANONICAL",
    "RCM_FAMILIES_CANONICAL",
    "ULLMANN_SPECIFIC_CANONICAL",
    "BUCHWALD_SPECIFIC_CANONICAL",
    "slugify_family",
    "canonical_family_label",
    "resolve_reaction_family",
    "apply_catalyst_override",
]

