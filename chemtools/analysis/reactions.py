"""
Reaction-level helpers that interface with the taxonomy registry.

This module centralises the family alias logic that was previously embedded in
``chemtools.router`` so that other components can consume the canonical
reaction identifiers without depending on the routing heuristics.
"""

from __future__ import annotations

import re
from typing import Optional, Dict, Set

from ._registry import get_registry

FAMILY_ALIAS_OVERRIDES: Dict[str, str] = {
    # C-N Coupling reactions
    "C_N_Coupling": "cn_coupling",
    "Buchwald_CN": "buchwald_hartwig_c_n",
    "Buchwald-Hartwig": "buchwald_hartwig_c_n",
    "Ullmann_CN": "ullmann_cn",
    "Chan_Lam_CN": "chan_lam",
    "chan_lam_cn": "chan_lam",  # Legacy alias
    
    # C-O Coupling reactions
    "C_O_Coupling": "co_coupling",
    "Ullmann_CO": "ullmann_ether",
    
    # C-S Coupling reactions
    "C_S_Coupling": "cs_coupling",
    
    # C-C Coupling reactions
    "Suzuki_CC": "suzuki_miyaura",
    "Suzuki": "suzuki_miyaura",
    "Suzuki-Miyaura": "suzuki_miyaura",
    "Negishi": "negishi",
    "Sonogashira_CC": "sonogashira",
    "Sonogashira": "sonogashira",
    "Stille": "stille",
    "Heck": "heck",
    
    # Amide Coupling
    "Amide_Coupling": "amide_coupling",
}

CN_FAMILIES_CANONICAL: Set[str] = {
    "cn_coupling",
    "buchwald_hartwig_c_n",
    "ullmann_cn",
    "chan_lam",
}
CO_FAMILIES_CANONICAL: Set[str] = {"co_coupling", "ullmann_ether"}
CS_FAMILIES_CANONICAL: Set[str] = {"cs_coupling"}
AMIDE_FAMILIES_CANONICAL: Set[str] = {"amide_coupling"}
SONOGASHIRA_FAMILIES_CANONICAL: Set[str] = {"sonogashira"}
SUZUKI_FAMILIES_CANONICAL: Set[str] = {"suzuki_miyaura", "suzuki_miyaura_in_situ"}
NEGISHI_FAMILIES_CANONICAL: Set[str] = {"negishi", "negishi_in_situ"}
STILLE_FAMILIES_CANONICAL: Set[str] = {"stille"}
HECK_FAMILIES_CANONICAL: Set[str] = {"heck"}
ULLMANN_SPECIFIC_CANONICAL: Set[str] = {"ullmann_cn"}
BUCHWALD_SPECIFIC_CANONICAL: Set[str] = {"buchwald_hartwig_c_n"}


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

    registry = get_registry()
    if registry:
        if family in registry.reaction_types:
            return family
        alias = registry.resolve_alias(family)
        if alias and alias.entity_type == "reaction_type":
            return alias.entity_id
        alias = registry.resolve_alias(family.lower())
        if alias and alias.entity_type == "reaction_type":
            return alias.entity_id
        slug = slugify_family(family)
        alias = registry.resolve_alias(slug)
        if alias and alias.entity_type == "reaction_type":
            return alias.entity_id
    return FAMILY_ALIAS_OVERRIDES.get(family)


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

    # Pd present suggests Buchwald-Hartwig unless strongly Ullmann-specific.
    if "Pd" in metals:
        if canonical in ULLMANN_SPECIFIC_CANONICAL and "Cu" not in metals:
            return canonical
        if canonical in CN_FAMILIES_CANONICAL or is_cn_coupling:
            return "buchwald_hartwig_c_n"
    # Cu without Pd leans toward Ullmann.
    if "Cu" in metals and canonical in CN_FAMILIES_CANONICAL:
        if canonical not in BUCHWALD_SPECIFIC_CANONICAL:
            return "ullmann_cn"
    return canonical


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
    "ULLMANN_SPECIFIC_CANONICAL",
    "BUCHWALD_SPECIFIC_CANONICAL",
    "slugify_family",
    "canonical_family_label",
    "resolve_reaction_family",
    "apply_catalyst_override",
]

