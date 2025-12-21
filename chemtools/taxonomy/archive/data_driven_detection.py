"""
Data-driven reaction type detection using taxonomy definitions.

This module provides reaction type detection based on the `reactants` field
in reaction_types.json, making the system extensible via JSON configuration
rather than hardcoded Python rules.

The detection works by:
1. Loading reaction types from reaction_types.json
2. Expanding `reactants` categories (e.g., "ArX*") to feature tokens
3. Matching detected features against required features

This replaces the hardcoded if/elif chain in _rule_based_detection().

Public Functions:
    detect_by_reactants(features: Dict[str, bool], catalysts: Set[str]) -> List[Match]
    load_reaction_requirements() -> Dict[str, ReactionRequirement]

Example:
    >>> features = {"ArBr_present": True, "boron_present": True}
    >>> catalysts = {"Pd"}
    >>> matches = detect_by_reactants(features, catalysts)
    >>> matches[0]
    Match(reaction_type="suzuki_miyaura", confidence=0.9, matched_features=...)
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple
import logging
from pathlib import Path
from functools import lru_cache
import json

from .reactant_mapper import (
    expand_category_to_features,
    expand_reactants_to_features,
    check_features_satisfy_reactants,
)

logger = logging.getLogger(__name__)


@dataclass
class ReactionRequirement:
    """Defines the reactant requirements for a reaction type."""
    id: str
    name: str
    category: Optional[str]
    reactants: List[str]  # Category IDs like ["ArX*", "ArB*"]
    catalysts: List[str]  # Required catalysts like ["Pd"]
    expanded_features: Dict[str, List[str]] = field(default_factory=dict)  # Cached expansion
    
    @property
    def requires_catalyst(self) -> bool:
        """Returns True if this reaction requires a metal catalyst."""
        metal_catalysts = set(self.catalysts) & _KNOWN_METAL_CATALYSTS
        return bool(metal_catalysts)


@dataclass
class DetectionMatch:
    """Represents a detection match result."""
    reaction_type: str
    name: str
    confidence: float
    matched_categories: Dict[str, List[str]]  # category -> matched features
    catalyst_match: bool
    reasons: List[str] = field(default_factory=list)


# ============================================================================
# Module-level caches
# ============================================================================
_REACTION_REQUIREMENTS_CACHE: Optional[Dict[str, ReactionRequirement]] = None

# Catalyst tokens in `reaction_types.json` are not uniformly structured; treat only
# known metal symbols as strict catalyst requirements/evidence.
_KNOWN_METAL_CATALYSTS: Set[str] = {"Pd", "Cu", "Ni", "Co"}


def _get_taxonomy_path() -> Path:
    """Get path to reaction_types.json."""
    return Path(__file__).parent / "data" / "reaction_types.json"


def load_reaction_requirements() -> Dict[str, ReactionRequirement]:
    """
    Load reaction type requirements from reaction_types.json.
    
    Returns:
        Dict mapping reaction_type_id to ReactionRequirement
    """
    global _REACTION_REQUIREMENTS_CACHE
    
    if _REACTION_REQUIREMENTS_CACHE is not None:
        return _REACTION_REQUIREMENTS_CACHE
    
    path = _get_taxonomy_path()
    
    if not path.exists():
        logger.warning(f"reaction_types.json not found at {path}")
        _REACTION_REQUIREMENTS_CACHE = {}
        return _REACTION_REQUIREMENTS_CACHE
    
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        _REACTION_REQUIREMENTS_CACHE = {}
        for entry in data:
            rxn_id = entry.get("id")
            if not rxn_id:
                continue
            
            reactants = entry.get("reactants", [])
            catalysts = entry.get("catalysts", [])
            
            req = ReactionRequirement(
                id=rxn_id,
                name=entry.get("name", rxn_id),
                category=entry.get("category"),
                reactants=reactants,
                catalysts=catalysts,
            )
            
            # Pre-expand reactants to features for faster matching
            if reactants:
                req.expanded_features = expand_reactants_to_features(reactants)
            
            _REACTION_REQUIREMENTS_CACHE[rxn_id] = req
        
        logger.debug(f"Loaded {len(_REACTION_REQUIREMENTS_CACHE)} reaction requirements")
        return _REACTION_REQUIREMENTS_CACHE
        
    except Exception as e:
        logger.warning(f"Failed to load reaction_types.json: {e}")
        _REACTION_REQUIREMENTS_CACHE = {}
        return _REACTION_REQUIREMENTS_CACHE


def _normalize_feature_keys(features: Dict[str, bool]) -> Dict[str, bool]:
    """
    Normalize feature keys to handle case variations.
    
    The calculable_features.json uses "ArBr_present" format, but
    _rule_hits may return different formats. This normalizes them.
    
    Also maps legacy keys like "aryl_halide" to category tokens.
    """
    # Create a case-insensitive mapping
    normalized = {}
    for key, value in features.items():
        # Keep original key
        normalized[key] = value
        # Also add lowercase version
        normalized[key.lower()] = value
    
    return normalized


def _map_legacy_features_to_tokens(features: Dict[str, bool]) -> Dict[str, bool]:
    """
    Map legacy feature names to reactant type tokens.
    
    The current _rule_hits returns keys like:
    - "aryl_halide" → should also check ArBr_present, ArCl_present, etc.
    - "boron" → should also check ArB(OH)2_present, ArB(OR)2_present, etc.
    - "terminal_alkyne" → should check terminal-alkyne_present
    
    This function expands legacy features to their modern token equivalents.
    """
    expanded = dict(features)  # Copy original
    
    # Legacy to category mappings
    legacy_to_category = {
        "aryl_halide": "ArX*",
        "vinyl_halide": "Vinyl-X",
        "boron": "ArB*",  # Also includes RB*
        "terminal_alkyne": "Alkyne",
        "grignard": "R-M",  # Grignard is part of organometallics
        "organozinc": "R-M",
        "organolithium": "R-M",
        "organostannane": "R-M",
        "nucleophile_n": "Amine",  # Primary/secondary amines
        "nucleophile_o": "Alcohol",
        "nucleophile_s": "Thiol",
        "alkene": "Alkene",
        "acid": "Carboxylic-acid",
        "alkyl_halide": "Alkyl-X",
    }
    
    # For each legacy feature that's True, set its category token
    # This helps with matching against expanded_features
    for legacy, category in legacy_to_category.items():
        if features.get(legacy):
            # Also set a category marker
            expanded[f"_legacy_{category}"] = True
    
    return expanded


def detect_by_reactants(
    features: Dict[str, bool],
    catalysts: Optional[Set[str]] = None,
    top_k: int = 5
) -> List[DetectionMatch]:
    """
    Detect reaction types based on detected features and catalysts.
    
    This is the main data-driven detection function that replaces
    the hardcoded if/elif chain.
    
    Args:
        features: Dict of detected features (from _rule_hits or feature detector)
        catalysts: Set of detected catalyst metals (e.g., {"Pd", "Cu"})
        top_k: Maximum number of matches to return
        
    Returns:
        List of DetectionMatch objects, sorted by confidence (highest first)
    """
    if catalysts is None:
        catalysts = set()
    
    requirements = load_reaction_requirements()
    
    if not requirements:
        return []
    
    # Normalize and expand features
    features = _normalize_feature_keys(features)
    expanded_features = _map_legacy_features_to_tokens(features)
    
    matches: List[DetectionMatch] = []
    
    for rxn_id, req in requirements.items():
        if not req.reactants:
            # Skip reactions without reactant definitions
            continue
        
        # Check if detected features satisfy required reactants
        satisfied, matched_cats = check_features_satisfy_reactants(
            expanded_features, 
            req.reactants
        )
        
        # Also check legacy feature mapping
        legacy_matches = _match_legacy_features(features, req.reactants)
        legacy_satisfied = all(bool(legacy_matches.get(cat_id)) for cat_id in req.reactants)
        
        if not satisfied and not legacy_satisfied:
            continue

        # Merge legacy matches into the main matched_cats so confidence and reasons
        # remain meaningful when the upstream feature detector provides legacy tokens.
        merged_matches = dict(matched_cats)
        for cat_id in req.reactants:
            if merged_matches.get(cat_id):
                continue
            legacy_hits = legacy_matches.get(cat_id) or []
            if legacy_hits:
                merged_matches[cat_id] = legacy_hits
        
        # Calculate confidence based on number of matched categories and catalysts
        base_confidence = 0.7
        
        # Boost for each matched category
        matched_count = sum(1 for hits in merged_matches.values() if hits)
        required_count = len(req.reactants)
        
        if required_count > 0:
            match_ratio = matched_count / required_count
            confidence = base_confidence + (0.2 * match_ratio)
        else:
            confidence = base_confidence
        
        # Check catalyst match
        catalyst_match = False
        required_metals = set(req.catalysts) & _KNOWN_METAL_CATALYSTS
        if catalysts and required_metals:
            if catalysts & required_metals:
                catalyst_match = True
                confidence = min(confidence + 0.1, 0.98)
            else:
                # Contradiction: detected metal catalyst doesn't match this reaction's typical metal.
                confidence = max(confidence - 0.1, 0.4)
        
        # Build reasons
        reasons = []
        for cat_id, hits in merged_matches.items():
            if hits:
                reasons.append(f"{cat_id}: {', '.join(hits[:3])}")
        if catalyst_match:
            reasons.append(f"catalyst: {catalysts & required_metals}")
        
        matches.append(DetectionMatch(
            reaction_type=rxn_id,
            name=req.name,
            confidence=confidence,
            matched_categories=merged_matches,
            catalyst_match=catalyst_match,
            reasons=reasons,
        ))
    
    # Sort by confidence (descending)
    matches.sort(key=lambda m: m.confidence, reverse=True)
    
    return matches[:top_k]


def _check_legacy_features(
    features: Dict[str, bool], 
    required_categories: List[str]
) -> bool:
    """
    Check if legacy feature names satisfy required categories.
    
    Maps legacy features like "aryl_halide" to categories like "ArX*".
    """
    matches = _match_legacy_features(features, required_categories)
    return all(bool(matches.get(cat_id)) for cat_id in required_categories)


def _match_legacy_features(
    features: Dict[str, bool],
    required_categories: List[str],
) -> Dict[str, List[str]]:
    """
    Return the legacy feature keys that satisfy each required category.

    This is used when the upstream feature detector emits coarse tokens (e.g.
    `aryl_halide`, `nucleophile_n`) instead of the fine-grained member tokens
    used by `reactant_mapper.expand_category_to_features(...)`.
    """
    category_to_legacy = {
        "ArX*": ["aryl_halide", "vinyl_halide", "triflate"],
        "HetAr-X": ["aryl_halide"],
        "Vinyl-X": ["vinyl_halide"],
        "ArB*": ["boron"],
        "RB*": ["boron"],
        "Alkene": ["alkene"],
        "Alkyne": ["terminal_alkyne"],
        "R-M": ["grignard", "organozinc", "organolithium", "organostannane"],
        "Amine": ["nucleophile_n"],
        "ArNH*": ["nucleophile_n"],
        "RNH2/R2NH": ["nucleophile_n"],
        "Alcohol": ["nucleophile_o", "alcohol"],
        "ROH": ["nucleophile_o", "alcohol"],
        "Thiol": ["nucleophile_s"],
        "RSH": ["nucleophile_s"],
        "Carboxylic-acid": ["acid"],
        "Alkyl-X": ["alkyl_halide"],
        "Aldehyde": ["carbonyl", "aldehyde"],
        "Ketone": ["carbonyl", "ketone"],
        "Acyl-electrophile": ["acyl_halide"],
    }

    matched: Dict[str, List[str]] = {}
    for cat_id in required_categories:
        legacy_names = category_to_legacy.get(cat_id, [])
        hits = [name for name in legacy_names if features.get(name)]
        matched[cat_id] = hits
    return matched


def get_reaction_requirement(reaction_type: str) -> Optional[ReactionRequirement]:
    """
    Get the requirement definition for a specific reaction type.
    
    Args:
        reaction_type: Reaction type ID (e.g., "suzuki_miyaura")
        
    Returns:
        ReactionRequirement or None if not found
    """
    requirements = load_reaction_requirements()
    return requirements.get(reaction_type)


def get_required_features_for_reaction(reaction_type: str) -> Dict[str, List[str]]:
    """
    Get the expanded feature requirements for a reaction type.
    
    Args:
        reaction_type: Reaction type ID (e.g., "suzuki_miyaura")
        
    Returns:
        Dict mapping category to list of feature tokens
        
    Example:
        >>> get_required_features_for_reaction("suzuki_miyaura")
        {
            "ArX*": ["ArBr_present", "ArCl_present", ...],
            "ArB*": ["ArB(OH)2_present", ...]
        }
    """
    req = get_reaction_requirement(reaction_type)
    if req:
        return req.expanded_features
    return {}


# ============================================================================
# Cache Management
# ============================================================================

def clear_caches():
    """Clear all module-level caches. Useful for testing or hot-reload."""
    global _REACTION_REQUIREMENTS_CACHE
    _REACTION_REQUIREMENTS_CACHE = None
    
    # Also clear reactant_mapper caches
    from . import reactant_mapper
    reactant_mapper.clear_caches()


if __name__ == "__main__":
    # Quick test
    print("=== Data-Driven Detection Test ===")
    
    # Load requirements
    reqs = load_reaction_requirements()
    print(f"\nLoaded {len(reqs)} reaction types")
    
    # Show a few examples
    for rxn_id in ["suzuki_miyaura", "sonogashira", "buchwald_hartwig_c_n", "heck"]:
        req = reqs.get(rxn_id)
        if req:
            print(f"\n{rxn_id}:")
            print(f"  Reactants: {req.reactants}")
            print(f"  Catalysts: {req.catalysts}")
            print(f"  Expanded features: {list(req.expanded_features.keys())}")
    
    # Test detection
    print("\n=== Detection Test ===")
    features = {
        "aryl_halide": True,
        "boron": True,
        "ArBr_present": True,
    }
    catalysts = {"Pd"}
    
    matches = detect_by_reactants(features, catalysts)
    print(f"\nWith features={features}, catalysts={catalysts}:")
    for m in matches[:3]:
        print(f"  {m.reaction_type}: {m.confidence:.2f} - {m.reasons}")
