"""
Centralized taxonomy mapping for all detection methods.

This module provides functions to map detection results (from SMARTS rules,
ML predictions, or heuristics) to canonical taxonomy IDs. All detection logic
must use these functions to ensure taxonomy alignment.
"""

from typing import Optional, Set, Dict
import re
import warnings

from .taxonomy import load_registry
from .taxonomy.registry import TaxonomyRegistry
from .analysis.reactions import canonical_family_label


_registry_cache: Optional[TaxonomyRegistry] = None
_registry_load_attempted = False


def get_taxonomy_registry() -> Optional[TaxonomyRegistry]:
    """
    Get cached taxonomy registry.
    
    Returns None if taxonomy cannot be loaded (allows graceful degradation).
    """
    global _registry_cache, _registry_load_attempted
    if not _registry_load_attempted:
        _registry_load_attempted = True
        try:
            _registry_cache = load_registry()
        except Exception as e:
            # Taxonomy validation errors should not break detection
            # Fall back to FAMILY_ALIAS_OVERRIDES in canonical_family_label
            warnings.warn(
                f"Failed to load taxonomy registry: {e}\n"
                "Detection will use fallback hardcoded mappings.",
                UserWarning,
                stacklevel=2
            )
    return _registry_cache


def resolve_to_taxonomy(
    raw_prediction: str,
    *,
    catalysts: Optional[Set[str]] = None,
    is_cn_coupling: bool = False,
    functional_groups: Optional[Dict[str, bool]] = None
) -> Optional[str]:
    """
    Resolve ANY detection result to canonical taxonomy ID.
    
    This is the ONLY function that should be used to convert
    detection results to taxonomy IDs. It provides robust mapping
    with three-layer strategy:
    
    1. Exact alias resolution via TaxonomyRegistry
    2. Keyword pattern matching for common variations
    3. Context-aware inference using catalysts + functional groups
    
    Args:
        raw_prediction: Raw string from SMARTS, ML, or heuristic
        catalysts: Set of detected catalyst metals (Pd, Cu, Ni, Co)
        is_cn_coupling: Whether C-N coupling signature detected
        functional_groups: Functional group hits from SMARTS detection
        
    Returns:
        Canonical taxonomy ID from reaction_types.json or None
        
    Examples:
        >>> resolve_to_taxonomy("Suzuki")
        "suzuki_miyaura"
        
        >>> resolve_to_taxonomy("C-N Coupling", catalysts={"Pd"})
        "buchwald_hartwig_c_n"
        
        >>> resolve_to_taxonomy("Suzuki coupling with boronic acids")
        "suzuki_miyaura"
        
        >>> resolve_to_taxonomy("Cross-coupling reaction", catalysts={"Pd"}, 
        ...                     functional_groups={"aryl_halide": True, "boron": True})
        "suzuki_miyaura"
    """
    if not raw_prediction:
        return None
    
    raw_prediction = raw_prediction.strip()
    if not raw_prediction:
        return None
    
    # Layer 1: Try direct taxonomy resolution (exact alias match)
    canonical = canonical_family_label(raw_prediction)
    if canonical:
        return canonical
    
    # If taxonomy unavailable, use fallback hardcoded mapping
    registry = get_taxonomy_registry()
    if not registry:
        # Fallback to basic mapping when taxonomy is unavailable
        fallback_map = {
            "suzuki_miyaura": "suzuki_miyaura",
            "suzuki": "suzuki_miyaura",
            "negishi": "negishi",
            "sonogashira": "sonogashira",
            "heck": "heck",
            "stille": "stille",
            "kumada": "kumada",
            "cn_coupling": "cn_coupling",
            "buchwald_hartwig_c_n": "buchwald_hartwig_c_n",
            "ullmann_cn": "ullmann_cn",
            "chan_lam": "chan_lam",
            "co_coupling": "co_coupling",
            "ullmann_ether": "ullmann_ether",
            "cs_coupling": "cs_coupling",
            "amide_coupling": "amide_coupling",
            "snar": "snar",
            "s_nar": "snar",
            "aromatic_nucleophilic_substitution": "snar",
            "reductive_amination": "reductive_amination",
            "esterification": "esterification",
            "grignard_addition": "grignard_addition",
            "organozinc_addition": "organozinc_addition",
            "organolithium_addition": "organolithium_addition",
            "hydrogenation": "hydrogenation",
            "carbonyl_reduction": "carbonyl_reduction",
            "alcohol_oxidation": "alcohol_oxidation",
            "epoxidation": "epoxidation",
            "e2_elimination": "e2_elimination",
            "williamson_ether": "williamson_ether",
            "finkelstein": "finkelstein",
            "nitrile_formation": "nitrile_formation",
            "hydroboration": "hydroboration",
            "diels_alder": "diels_alder",
            "michael_addition": "michael_addition",
            "claisen_condensation": "claisen_condensation",
            "aldol_condensation": "aldol_condensation",
            "wittig": "wittig",
        }
        if raw_prediction.lower() in fallback_map:
            return fallback_map[raw_prediction.lower()]
        # Try without underscores
        cleaned = raw_prediction.lower().replace("-", "_").replace(" ", "_")
        if cleaned in fallback_map:
            return fallback_map[cleaned]
    
    # Layer 2: Try case-insensitive and slug variations
    lower_pred = raw_prediction.lower()
    canonical = canonical_family_label(lower_pred)
    if canonical:
        return canonical
    
    # Try slugified version
    slug = re.sub(r"[^0-9a-z]+", "_", lower_pred).strip("_")
    canonical = canonical_family_label(slug)
    if canonical:
        return canonical
    
    # Layer 3: Keyword pattern matching for common ML variations
    # (rxn-insight returns unpredictable names like "Suzuki coupling" vs "Cross-coupling reaction")
    keyword_map = {
        "suzuki": "suzuki_miyaura",
        "negishi": "negishi",
        "stille": "stille",
        "heck": "heck",
        "sonogashira": "sonogashira",
        "buchwald": "buchwald_hartwig_c_n",
        "hartwig": "buchwald_hartwig_c_n",
        "ullmann": "ullmann_cn",
        "chan-lam": "chan_lam",
        "chan lam": "chan_lam",
        "grignard": "grignard_addition",
        "kumada": "kumada",
        "amide": "amide_coupling",
        "snar": "snar",
        "s_nar": "snar",
        "aromatic nucleophilic": "snar",
        "reductive amination": "reductive_amination",
        "esterification": "esterification",
        "hydrogenation": "hydrogenation",
        "diels-alder": "diels_alder",
        "diels alder": "diels_alder",
        "michael": "michael_addition",
        "aldol": "aldol_condensation",
        "claisen": "claisen_condensation",
        "wittig": "wittig",
        "reduction": "carbonyl_reduction",
        "oxidation": "alcohol_oxidation",
    }
    
    for keyword, family_id in keyword_map.items():
        if keyword in lower_pred:
            # Verify it exists in taxonomy (if available)
            registry = get_taxonomy_registry()
            if registry and registry.get_reaction_type(family_id):
                return family_id
            elif not registry:
                # Taxonomy unavailable, trust the mapping
                return family_id
    
    # Layer 4: Context-aware inference for ambiguous predictions
    # Use catalysts + functional groups to disambiguate
    catalysts = catalysts or set()
    functional_groups = functional_groups or {}
    
    # Handle generic "C-C coupling" or "cross-coupling" predictions
    if any(term in lower_pred for term in ["c-c coupling", "c-c bond", "cross-coupling", "cross coupling"]):
        # Use functional groups to determine specific type
        if functional_groups.get("boron"):
            return "suzuki_miyaura"
        elif functional_groups.get("terminal_alkyne"):
            return "sonogashira"
        elif functional_groups.get("grignard"):
            return "kumada"
        elif functional_groups.get("organozinc"):
            return "negishi"
        elif functional_groups.get("alkene"):
            return "heck"
        # Generic C-C coupling without more context
        return None
    
    # Handle generic "C-N coupling" or "amination" predictions
    if (is_cn_coupling or 
        any(term in lower_pred for term in ["c-n coupling", "c-n bond", "amination", "n-arylation"])):
        # Use catalyst to determine specific type
        if "Pd" in catalysts:
            return "buchwald_hartwig_c_n"
        elif "Cu" in catalysts:
            return "ullmann_cn"
        # Generic C-N coupling without catalyst info
        return "cn_coupling"
    
    # Handle generic "C-O coupling" predictions
    if any(term in lower_pred for term in ["c-o coupling", "c-o bond", "etherification"]):
        if "Cu" in catalysts:
            return "ullmann_ether"
        return "co_coupling"
    
    # Handle generic "C-S coupling" predictions
    if any(term in lower_pred for term in ["c-s coupling", "c-s bond", "thioetherification"]):
        return "cs_coupling"
    
    # Conservative fallback: return None if uncertain
    # Better to return None than misclassify
    return None


def calculate_confidence_adjustment(
    raw_prediction: str,
    mapped_family: Optional[str],
    base_confidence: float
) -> float:
    """
    Adjust confidence based on mapping certainty.
    
    Exact alias matches get no penalty, keyword matches get small penalty,
    context-based inference gets larger penalty.
    
    Args:
        raw_prediction: Original ML/rule prediction
        mapped_family: Mapped canonical ID (or None)
        base_confidence: Original confidence from detector
        
    Returns:
        Adjusted confidence score (0.0-1.0)
    """
    if not mapped_family:
        return 0.0
    
    # Check if exact alias match
    canonical = canonical_family_label(raw_prediction)
    if canonical == mapped_family:
        # Exact match - no penalty
        return base_confidence * 0.98  # Tiny penalty for being mapped
    
    # Check if keyword match
    lower_pred = raw_prediction.lower()
    if mapped_family.replace("_", " ") in lower_pred or mapped_family.replace("_", "-") in lower_pred:
        # Keyword match - small penalty
        return base_confidence * 0.90
    
    # Check if contains family name
    family_keywords = {
        "suzuki": ["suzuki", "miyaura"],
        "buchwald": ["buchwald", "hartwig"],
        "ullmann": ["ullmann", "goldberg"],
        "sonogashira": ["sonogashira"],
        "negishi": ["negishi"],
        "stille": ["stille"],
        "heck": ["heck"],
    }
    
    family_base = mapped_family.split("_")[0]
    keywords = family_keywords.get(family_base, [family_base])
    if any(kw in lower_pred for kw in keywords):
        # Partial match - medium penalty
        return base_confidence * 0.85
    
    # Context-based inference - larger penalty
    return base_confidence * 0.70


def get_mapping_method(raw_prediction: str, mapped_family: Optional[str]) -> str:
    """
    Track HOW we mapped the raw prediction to taxonomy.
    
    Useful for debugging and improving taxonomy coverage.
    
    Returns:
        One of: "exact_alias", "keyword_match", "context_inference", "failed"
    """
    if not mapped_family:
        return "failed"
    
    # Check exact alias
    canonical = canonical_family_label(raw_prediction)
    if canonical == mapped_family:
        return "exact_alias"
    
    # Check keyword match
    lower_pred = raw_prediction.lower()
    if mapped_family.replace("_", " ") in lower_pred or mapped_family.replace("_", "-") in lower_pred:
        return "keyword_match"
    
    # Must be context inference
    return "context_inference"
