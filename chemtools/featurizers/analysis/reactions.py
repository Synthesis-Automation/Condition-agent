"""
Reaction-level helpers that interface with the taxonomy registry.

This module centralises the family alias logic that was previously embedded in
``chemtools.router`` so that other components can consume the canonical
reaction identifiers without depending on the routing heuristics.

IMPORTANT: Reaction types should be defined in taxonomy/data/reaction_types.v4.0.json.
This alias mapping is for legacy compatibility only. New reactions should be added
to the taxonomy JSON, not to FAMILY_ALIAS_OVERRIDES.

The system is taxonomy-driven: reaction detection uses the definitions in reaction_types
JSON, and this module only provides backward-compatible alias resolution.
"""

from __future__ import annotations

import re
from typing import Optional, Dict, Set

from ...taxonomy import reaction_catalog as _reaction_catalog

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
    "ChanLam_Narylation": "c_n_cross_coupling",
    "UllmannGoldberg_Narylation": "c_n_cross_coupling",
    
    # C-C Coupling reactions
    "Suzuki_CC": "suzuki_miyaura",
    "Suzuki": "suzuki_miyaura",
    "Suzuki-Miyaura": "suzuki_miyaura",
    "Sonogashira_CC": "sonogashira",
    "Sonogashira_coupling": "sonogashira",
    "Heck": "heck",
    "HeckMizoroki_coupling": "heck",
    # Note: Stille, Negishi, Kumada, Hiyama are distinct reactions defined in taxonomy
    # They should NOT map to suzuki_miyaura - each uses different organometallic reagents
    "Stille_coupling": "Stille",
    "Negishi_coupling": "Negishi",
    "Kumada_coupling": "Kumada",
    "Hiyama_coupling": "Hiyama",
    "Miyaura_borylation": "suzuki_miyaura",
    "Olefin_metathesis": "ring_closing_metathesis",
    "RCM": "ring_closing_metathesis",
    
    # Amide Coupling
    "Amide_Coupling": "amide_coupling",
    "Amide_formation": "amide_coupling",
    "amide_formation": "amide_coupling",
    "CDI-mediated_amidation": "amide_coupling",
    "Weinreb_amide": "amide_coupling",
    
    # SNAr reactions (legacy -> v2)
    "S_NAr": "snar_cn",
    "SNAr": "snar_cn",
    "snar": "snar_cn",
    "s_nar": "snar_cn",
    "Aromatic_Nucleophilic_Substitution": "snar_cn",
    "aromatic_nucleophilic_substitution": "snar_cn",
    "SNAr_amination": "snar_cn",
    "SNAr-CN": "snar_cn",
    "SNAr-CO": "snar_cn",
    "SNAr-CS": "snar_cn",
    
    # Esterification
    "Esterification": "esterification",
    "Ester_formation": "esterification",
    "Fischer_esterification": "esterification",
    "Steglich_esterification": "esterification",

    # Reductive Amination
    "Reductive_amination": "reductive_amination",

    # C-O Coupling
    "C_O_Coupling": "c_o_coupling",

    # Oxidation
    "Swern_oxidation": "oxidation_alcohol_to_carbonyl",
    "Jones_oxidation": "oxidation_alcohol_to_carbonyl",
    "DessMartin_periodinane_DMP_Alcohols__aldehydesketones": "oxidation_alcohol_to_carbonyl",
    "Oxidation-BAIB": "oxidation_alcohol_to_carbonyl",
    "Oxidation_IBX": "oxidation_alcohol_to_carbonyl",
    "Oxidation_of_Aromatic_Side_Chains": "oxidation_alcohol_to_carbonyl",
    "Oxidation_of_Arylmethanes_to_Aldehydes": "oxidation_alcohol_to_carbonyl",
    "Riley_oxidation": "oxidation_alcohol_to_carbonyl",
    "TEMPObleach_Oxidations_Primary_alcohols__aldehydesacids": "oxidation_alcohol_to_carbonyl",
    "TPAPNMO_Catalytic_Ru_oxidations": "oxidation_alcohol_to_carbonyl",
    "Wacker_oxidation_Alkene__ketone_PdCl2_CuCl_O2": "oxidation_alcohol_to_carbonyl",
    "Benzylic_oxidation": "oxidation_alcohol_to_carbonyl",

    # Reduction
    "NaBH4_carbonyl_reductions": "reduction_carbonyl_to_alcohol",
    "Reduction_LAH_LiAlH4": "reduction_carbonyl_to_alcohol",
    "Reduction_NaBHOAc3_NaBH3CN": "reduction_carbonyl_to_alcohol",
    "Reduction-Borane": "reduction_carbonyl_to_alcohol",
    "DIBALH_Partial_reductions": "reduction_carbonyl_to_alcohol",
    "Clemmensen_reduction": "reduction_carbonyl_to_alcohol",
    "Rosenmund_reduction": "reduction_carbonyl_to_alcohol",
    "WolffKishner": "reduction_carbonyl_to_alcohol",
    "Catalytic_hydrogenation": "reduction_nitro_to_amine",
    "Staudinger_reduction": "reduction_nitro_to_amine",
    "Birch_reduction": "reduction_carbonyl_to_alcohol",

    # Halogenation
    "NBS_NCS_NIS_halogenation": "halogenation_aromatic",
    "Sandmeyer_reactions": "halogenation_aromatic",
    "Appel_halogenation": "halogenation_aromatic",
    "Aromatic_Halogen_Exchange": "halogenation_aromatic",
    "BalzSchiemann": "halogenation_aromatic",
    "Chlorination_SOCl2_oxalyl_chloride": "halogenation_aromatic",
    "Deoxy_fluorination": "halogenation_aromatic",
    "Electrophilic_fluorination": "halogenation_aromatic",
    "Trifluoromethylation": "halogenation_aromatic",
    "Difluoromethylation": "halogenation_aromatic",
    "Aliphatic_Halide_Exchange": "halogenation_aromatic",
}

CN_FAMILIES_CANONICAL: Set[str] = {"c_n_cross_coupling", "snar_cn"}
CO_FAMILIES_CANONICAL: Set[str] = {"c_o_coupling"}
CS_FAMILIES_CANONICAL: Set[str] = {"c_s_coupling"}
AMIDE_FAMILIES_CANONICAL: Set[str] = {"amide_coupling"}
SONOGASHIRA_FAMILIES_CANONICAL: Set[str] = {"sonogashira", "Sonogashira"}
SUZUKI_FAMILIES_CANONICAL: Set[str] = {"suzuki_miyaura", "Suzuki_miyaura"}
NEGISHI_FAMILIES_CANONICAL: Set[str] = {"Negishi"}
STILLE_FAMILIES_CANONICAL: Set[str] = {"Stille"}
KUMADA_FAMILIES_CANONICAL: Set[str] = {"Kumada"}
HIYAMA_FAMILIES_CANONICAL: Set[str] = {"Hiyama"}
HECK_FAMILIES_CANONICAL: Set[str] = {"heck", "Heck"}
RCM_FAMILIES_CANONICAL: Set[str] = {"ring_closing_metathesis"}
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

