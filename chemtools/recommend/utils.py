"""
Utility functions for reaction recommendation system.

Helper functions for family normalization, reactant analysis, and constraint handling.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

from .. import constraints
from ..synthon import select_electrophile_nucleophile
from ..taxonomy import reaction_catalog as _reaction_catalog


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
    "snar_cn": "C_N_Coupling",
    "snar_co": "C_O_Coupling",
    "snar_cs": "C_S_Coupling",
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
    "Sonogashira_CC": "Sonogashira_coupling",
    "sonogashira": "Sonogashira_coupling",
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
    "Heck": "HeckMizoroki_coupling",
    "amide_coupling": "Amide_formation",
}

_ALIAS_CASEFOLD = {k.lower(): v for k, v in _FAMILY_ALIASES.items()}

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


@lru_cache(maxsize=1)
def _dataset_family_index() -> Dict[str, str]:
    base = Path(__file__).resolve().parent.parent.parent / "data" / "reaction_dataset"
    if not base.exists():
        return {}
    return {path.stem.lower(): path.stem for path in base.glob("*.jsonl")}


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
    if not fam:
        return "Unknown"

    resolved = _reaction_catalog.resolve_reaction_type(fam)
    if resolved is None:
        resolved = _reaction_catalog.resolve_reaction_type(fam.lower())
    if resolved:
        fam = resolved

    mapped = _FAMILY_ALIASES.get(fam) or _ALIAS_CASEFOLD.get(fam.lower())
    candidate = mapped or fam

    dataset_index = _dataset_family_index()
    if dataset_index:
        direct = dataset_index.get(candidate.lower())
        if direct:
            return direct
        direct = dataset_index.get(fam.lower())
        if direct:
            return direct
        variants = {
            candidate.replace(" ", "_"),
            candidate.replace(" ", "-"),
            candidate.replace("_", "-"),
            candidate.replace("-", "_"),
        }
        for variant in variants:
            direct = dataset_index.get(variant.lower())
            if direct:
                return direct

    return candidate


def _taxonomy_family_label(family: str) -> Optional[str]:
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
    
    Uses taxonomy-driven synthon assignment with functional-group fallback.
    
    Args:
        reactants: List of reactant SMILES
        
    Returns:
        (electrophile_smiles, nucleophile_smiles)
    """
    return select_electrophile_nucleophile(reactants)


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
