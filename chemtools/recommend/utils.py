"""
Utility functions for reaction recommendation system.

Helper functions for family normalization, reactant analysis, and constraint handling.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, Any, List, Tuple, Optional

from .. import constraints
from ..taxonomy import load_registry
from ..taxonomy.v2 import reaction_catalog as _reaction_catalog


# Family name normalization
_FAMILY_ALIASES: Dict[str, str] = {
    # Suzuki variants -> canonical dataset name
    "Suzuki_CC": "Suzuki",
    "Suzuki_Coupling": "Suzuki",
    "Suzuki": "Suzuki",
    "suzuki_miyaura": "Suzuki",
    "negishi": "Negishi",
    "stille": "Stille",
    "heck": "Heck",
    # C-N Coupling variants -> unified canonical dataset name
    "Ullmann C?N": "C_N_Coupling",
    "Ullmann_CN": "C_N_Coupling",
    "Buchwald C?N": "C_N_Coupling",
    "Buchwald_CN": "C_N_Coupling",
    "C_N_Coupling_Cu": "C_N_Coupling",  # Legacy copper-specific
    "C_N_Coupling_Pd": "C_N_Coupling",  # Legacy palladium-specific
    "C_N_Coupling_Ni": "C_N_Coupling",  # Legacy nickel-specific
    "C_N_Coupling_Cu_Ullmann": "C_N_Coupling",
    "C_N_Coupling_Pd_Buchwald": "C_N_Coupling",
    "C_N_Coupling": "C_N_Coupling",
    "cn_coupling": "C_N_Coupling",
    "buchwald_hartwig_c_n": "C_N_Coupling",
    "ullmann_cn": "C_N_Coupling",
    "chan_lam": "C_N_Coupling",
    "chan_lam_cn": "C_N_Coupling",  # Legacy alias
    "c_n_cross_coupling": "C_N_Coupling",
    "snar_cn": "SNAr-CN",
    # C-O Coupling variants -> canonical dataset name
    "C_O_Coupling": "C_O_Coupling",
    "Ullmann_CO": "C_O_Coupling",
    "Ullmann C?O": "C_O_Coupling",
    "co_coupling": "C_O_Coupling",
    "ullmann_ether": "C_O_Coupling",
    # C-S Coupling variants -> canonical dataset name
    "C_S_Coupling": "C_S_Coupling",
    "Ullmann_CS": "C_S_Coupling",
    "Ullmann C?S": "C_S_Coupling",
    "cs_coupling": "C_S_Coupling",
    # Sonogashira variants
    "Sonogashira_CC": "Sonogashira_CC",
    "sonogashira": "Sonogashira_CC",
    "sonogashira_coupling": "Sonogashira_coupling",
    # Amide variants -> canonical dataset name
    "Amide_Coupling": "Amide_formation",
    "Amide_formation": "Amide_formation",
    "amide_coupling": "Amide_formation",
    "amide_formation": "Amide_formation",
    # v2 IDs mapped to dataset names
    "suzuki_miyaura": "Suzuki",
    "sonogashira": "Sonogashira_coupling",
    "heck": "HeckMizoroki_coupling",
    "amide_coupling": "Amide_formation",
}

_FAMILY_LABELS: Dict[str, str] = {
    "c_n_cross_coupling": "C–N Cross-Coupling",
    "snar_cn": "SNAr Amination",
    "sn2_substitution": "SN2 Substitution",
    "suzuki_miyaura": "Suzuki–Miyaura Coupling",
    "sonogashira": "Sonogashira Coupling",
    "heck": "Heck Reaction",
    "amide_coupling": "Amide Formation",
    "esterification": "Esterification",
    "acylation_acyl_halide_amide": "Acyl Halide Amidation",
    "acylation_acyl_halide_ester": "Acyl Halide Esterification",
    # Legacy aliases
    "cn_coupling": "C–N Coupling",
    "buchwald_hartwig_c_n": "Buchwald–Hartwig C–N Coupling",
    "ullmann_cn": "Ullmann C–N Coupling",
    "chan_lam": "Chan–Lam Coupling",
    "chan_lam_cn": "Chan–Lam Coupling",
    "co_coupling": "C–O Coupling",
    "ullmann_ether": "Ullmann Ether Synthesis",
    "cs_coupling": "C–S Coupling",
    "negishi": "Negishi Coupling",
    "stille": "Stille Coupling",
}


def canonical_family(family: str | None) -> str:
    """
    Normalize family name to canonical form.
    
    Args:
        family: Raw family name
        
    Returns:
        Canonical family name or "Unknown"
    """
    if not family:
        return "Unknown"
    fam = str(family).strip()
    return _FAMILY_ALIASES.get(fam, fam)


@lru_cache(maxsize=1)
def _cached_registry():
    try:
        return load_registry()
    except Exception:
        return None


def _taxonomy_family_label(family: str) -> Optional[str]:
    registry = _cached_registry()
    if registry is None:
        return _v2_family_label(family)
    record = registry.get_reaction_type(family)
    if record is not None:
        return record.name
    alias = registry.resolve_alias(family)
    if alias and alias.entity_type == "reaction_type":
        target = registry.get_reaction_type(alias.entity_id)
        if target is not None:
            return target.name
    return _v2_family_label(family)


def _v2_family_label(family: str) -> Optional[str]:
    resolved = _reaction_catalog.resolve_reaction_type(family)
    if resolved:
        record = _reaction_catalog.get_reaction_type(resolved)
        if record is not None:
            return record.name
    return None


def friendly_family_label(family: str | None) -> Optional[str]:
    """
    Return a human-friendly label for a canonical reaction family.
    """
    if not family:
        return None
    fam = str(family).strip()
    if not fam:
        return None

    # Try taxonomy lookups (exact, alias, and simple variants)
    label = _taxonomy_family_label(fam)
    if label:
        return label

    variants = {
        fam.lower(),
        fam.replace("_", " "),
        fam.replace("_", "-"),
        fam.replace("-", " "),
    }
    for candidate in variants:
        candidate = candidate.strip()
        if not candidate:
            continue
        label = _taxonomy_family_label(candidate)
        if label:
            return label

    canonical = canonical_family(fam)
    label = _FAMILY_LABELS.get(canonical) or _FAMILY_LABELS.get(fam)
    if label:
        return label

    # Fallback: prettify canonical ID.
    return canonical.replace("_", " ").replace("-", " ").title()


def pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    """
    Pick electrophile and nucleophile from reactant list.
    
    Uses simple heuristic: looks for halides or sulfonates.
    
    Args:
        reactants: List of reactant SMILES
        
    Returns:
        (electrophile_smiles, nucleophile_smiles)
    """
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return (
            ("br" in t) or ("cl" in t) or (" i" in t)
            or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
        )
    
    if not reactants:
        return "", ""
    if len(reactants) == 1:
        return reactants[0], ""
    
    r0, r1 = reactants[0], reactants[1]
    if is_electrophile(r0):
        return r0, r1
    if is_electrophile(r1):
        return r1, r0
    
    return r0, r1


def median(vals: List[float]) -> float | None:
    """
    Compute median of numeric list.
    
    Args:
        vals: List of numbers
        
    Returns:
        Median value or None if empty
    """
    xs = [float(v) for v in vals if isinstance(v, (int, float))]
    if not xs:
        return None
    xs.sort()
    n = len(xs)
    mid = n // 2
    if n % 2 == 1:
        return xs[mid]
    return 0.5 * (xs[mid - 1] + xs[mid])


def pick_with_constraints(
    cands: List[str],
    rules: Dict[str, Any]
) -> Tuple[str | None, Dict[str, Any]]:
    """
    Apply constraints to candidate list and pick first allowed item.
    
    Args:
        cands: List of candidate reagent names
        rules: Constraint rules dict
        
    Returns:
        (picked_item, filter_result_dict)
    """
    if not cands:
        return None, {"allowed": [], "blocked": []}
    if not rules:
        return cands[0], {"allowed": cands, "blocked": []}
    
    out = constraints.apply_filter(cands, rules)
    allowed = out.get("allowed") or []
    return (allowed[0] if allowed else None), out
