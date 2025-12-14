"""
Reactant category to feature token mapper.

This module provides utilities to expand reactant categories (e.g., "ArX*")
into their constituent feature tokens (e.g., ["ArBr_present", "ArCl_present", ...]).

The mapping is driven by:
- chemtools/taxonomy/data/reactant_types.json: Defines reactant categories and members
- chemtools/taxonomy/data/calculable_features.json: Maps members to feature tokens

Public Functions:
    load_reactant_taxonomy() -> Dict[str, List[str]]
    expand_category_to_features(category_id: str) -> List[str]
    expand_reactants_to_features(reactant_ids: List[str]) -> Dict[str, List[str]]
    check_features_satisfy_reactants(detected: Dict[str, bool], required: List[str]) -> Tuple[bool, Dict]

Example:
    >>> expand_category_to_features("ArX*")
    ["ArBr_present", "ArCl_present", "ArI_present", "ArOTf_present", ...]
    
    >>> check_features_satisfy_reactants(
    ...     detected={"ArBr_present": True, "boron_present": True},
    ...     required=["ArX*", "ArB*"]
    ... )
    (True, {"ArX*": ["ArBr_present"], "ArB*": []})  # boron_present belongs to ArB*
"""

from typing import Any, Dict, List, Optional, Tuple
import logging
from pathlib import Path
from functools import lru_cache
import json

logger = logging.getLogger(__name__)

# ============================================================================
# Module-level caches
# ============================================================================
_REACTANT_TYPES_CACHE: Optional[Dict[str, Dict[str, Any]]] = None
_MEMBER_TO_TOKEN_CACHE: Optional[Dict[str, str]] = None
_CATEGORY_TO_MEMBERS_CACHE: Optional[Dict[str, List[str]]] = None


def _get_data_path(filename: str) -> Path:
    """Get path to data file in taxonomy/data directory."""
    return Path(__file__).parent / "data" / filename


def _load_reactant_types() -> Dict[str, Dict[str, Any]]:
    """
    Load reactant_types.json and index by category ID.
    
    Returns:
        Dict mapping category_id -> category definition
        Example: {"ArX*": {"id": "ArX*", "name": "aryl halide", "members": [...]}, ...}
    """
    global _REACTANT_TYPES_CACHE
    
    if _REACTANT_TYPES_CACHE is not None:
        return _REACTANT_TYPES_CACHE
    
    path = _get_data_path("reactant_types.json")
    
    if not path.exists():
        logger.warning(f"reactant_types.json not found at {path}")
        _REACTANT_TYPES_CACHE = {}
        return _REACTANT_TYPES_CACHE
    
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        _REACTANT_TYPES_CACHE = {}
        for entry in data:
            cat_id = entry.get("id")
            if cat_id:
                _REACTANT_TYPES_CACHE[cat_id] = entry
        
        logger.debug(f"Loaded {len(_REACTANT_TYPES_CACHE)} reactant categories")
        return _REACTANT_TYPES_CACHE
        
    except Exception as e:
        logger.warning(f"Failed to load reactant_types.json: {e}")
        _REACTANT_TYPES_CACHE = {}
        return _REACTANT_TYPES_CACHE


def _build_member_to_token_map() -> Dict[str, str]:
    """
    Build mapping from reactant member ID to feature token.
    
    Scans calculable_features.json for entries with reactant_metadata.reactant_member
    and maps member ID -> token name.
    
    Returns:
        Dict: {"ArBr": "ArBr_present", "ArCl": "ArCl_present", ...}
    """
    global _MEMBER_TO_TOKEN_CACHE
    
    if _MEMBER_TO_TOKEN_CACHE is not None:
        return _MEMBER_TO_TOKEN_CACHE
    
    path = _get_data_path("calculable_features.json")
    
    if not path.exists():
        logger.warning(f"calculable_features.json not found at {path}")
        _MEMBER_TO_TOKEN_CACHE = {}
        return _MEMBER_TO_TOKEN_CACHE
    
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        features = data.get("features", []) if isinstance(data, dict) else data
        
        _MEMBER_TO_TOKEN_CACHE = {}
        for feat in features:
            meta = feat.get("reactant_metadata", {})
            member = meta.get("reactant_member")
            token = feat.get("token")
            if member and token:
                _MEMBER_TO_TOKEN_CACHE[member] = token
        
        logger.debug(f"Built member->token map with {len(_MEMBER_TO_TOKEN_CACHE)} entries")
        return _MEMBER_TO_TOKEN_CACHE
        
    except Exception as e:
        logger.warning(f"Failed to build member->token map: {e}")
        _MEMBER_TO_TOKEN_CACHE = {}
        return _MEMBER_TO_TOKEN_CACHE


def _build_category_to_members_map() -> Dict[str, List[str]]:
    """
    Build mapping from reactant category ID to list of member IDs.
    
    Returns:
        Dict: {"ArX*": ["ArBr", "ArCl", "ArI", ...], "ArB*": ["ArB(OH)2", ...], ...}
    """
    global _CATEGORY_TO_MEMBERS_CACHE
    
    if _CATEGORY_TO_MEMBERS_CACHE is not None:
        return _CATEGORY_TO_MEMBERS_CACHE
    
    reactant_types = _load_reactant_types()
    
    _CATEGORY_TO_MEMBERS_CACHE = {}
    for cat_id, cat_def in reactant_types.items():
        members = cat_def.get("members", [])
        member_ids = [m.get("id") for m in members if m.get("id")]
        _CATEGORY_TO_MEMBERS_CACHE[cat_id] = member_ids
    
    return _CATEGORY_TO_MEMBERS_CACHE


# ============================================================================
# Public API
# ============================================================================

def load_reactant_taxonomy() -> Dict[str, Dict[str, Any]]:
    """
    Load the full reactant taxonomy.
    
    Returns:
        Dict mapping category_id to category definition with members.
    """
    return _load_reactant_types()


def get_category_members(category_id: str) -> List[str]:
    """
    Get the member IDs for a reactant category.
    
    Args:
        category_id: Category ID like "ArX*" or "ArB*"
        
    Returns:
        List of member IDs: ["ArBr", "ArCl", ...]
    """
    cat_to_members = _build_category_to_members_map()
    return cat_to_members.get(category_id, [])


def expand_category_to_features(category_id: str) -> List[str]:
    """
    Expand a reactant category to its feature tokens.
    
    Args:
        category_id: Category ID like "ArX*" or "ArB*"
        
    Returns:
        List of feature tokens: ["ArBr_present", "ArCl_present", ...]
        Returns empty list if category not found.
        
    Example:
        >>> expand_category_to_features("ArX*")
        ["ArBr_present", "ArCl_present", "ArI_present", "ArF_present",
         "ArOSO2R_present", "ArOMs_present", "ArOTf_present", "ArOTs_present"]
    """
    member_ids = get_category_members(category_id)
    member_to_token = _build_member_to_token_map()
    
    tokens = []
    for member_id in member_ids:
        token = member_to_token.get(member_id)
        if token:
            tokens.append(token)
        else:
            # Fallback: try "<member_id>_present" pattern
            fallback = f"{member_id}_present"
            tokens.append(fallback)
    
    return tokens


def expand_reactants_to_features(reactant_ids: List[str]) -> Dict[str, List[str]]:
    """
    Expand multiple reactant categories to their feature tokens.
    
    Args:
        reactant_ids: List of category IDs like ["ArX*", "ArB*"]
        
    Returns:
        Dict mapping category_id to list of feature tokens
        
    Example:
        >>> expand_reactants_to_features(["ArX*", "ArB*"])
        {
            "ArX*": ["ArBr_present", "ArCl_present", ...],
            "ArB*": ["ArB(OH)2_present", "ArB(OR)2_present", ...]
        }
    """
    result = {}
    for cat_id in reactant_ids:
        result[cat_id] = expand_category_to_features(cat_id)
    return result


def check_features_satisfy_reactants(
    detected: Dict[str, bool], 
    required: List[str]
) -> Tuple[bool, Dict[str, List[str]]]:
    """
    Check if detected features satisfy the required reactant categories.
    
    For a reaction type to match, at least one feature from each required
    category must be detected.
    
    Args:
        detected: Dict of detected features, e.g., {"ArBr_present": True, "boron_present": True}
        required: List of required category IDs, e.g., ["ArX*", "ArB*"]
        
    Returns:
        Tuple of (satisfied: bool, matched_features: Dict[str, List[str]])
        - satisfied: True if all required categories have at least one detected feature
        - matched_features: Dict mapping category to list of matched feature tokens
        
    Example:
        >>> check_features_satisfy_reactants(
        ...     detected={"ArBr_present": True, "ArB(OH)2_present": True},
        ...     required=["ArX*", "ArB*"]
        ... )
        (True, {"ArX*": ["ArBr_present"], "ArB*": ["ArB(OH)2_present"]})
    """
    # Get only features that are True
    active_features = {k for k, v in detected.items() if v}
    
    matched = {}
    all_satisfied = True
    
    for cat_id in required:
        category_tokens = set(expand_category_to_features(cat_id))
        hits = list(active_features & category_tokens)
        matched[cat_id] = hits
        
        if not hits:
            all_satisfied = False
    
    return all_satisfied, matched


def get_all_categories() -> List[str]:
    """
    Get list of all available reactant category IDs.
    
    Returns:
        List of category IDs: ["ArX*", "ArB*", "Alkene", ...]
    """
    return list(_load_reactant_types().keys())


def get_category_info(category_id: str) -> Optional[Dict[str, Any]]:
    """
    Get full information about a reactant category.
    
    Args:
        category_id: Category ID like "ArX*"
        
    Returns:
        Category definition dict or None if not found.
    """
    return _load_reactant_types().get(category_id)


# ============================================================================
# Cache Management
# ============================================================================

def clear_caches():
    """Clear all module-level caches. Useful for testing or hot-reload."""
    global _REACTANT_TYPES_CACHE, _MEMBER_TO_TOKEN_CACHE, _CATEGORY_TO_MEMBERS_CACHE
    _REACTANT_TYPES_CACHE = None
    _MEMBER_TO_TOKEN_CACHE = None
    _CATEGORY_TO_MEMBERS_CACHE = None


if __name__ == "__main__":
    # Quick test
    print("=== Reactant Mapper Test ===")
    
    # Test category expansion
    for cat_id in ["ArX*", "ArB*", "Alkene", "Alkyne"]:
        tokens = expand_category_to_features(cat_id)
        print(f"\n{cat_id} -> {len(tokens)} tokens:")
        for t in tokens[:5]:
            print(f"  - {t}")
        if len(tokens) > 5:
            print(f"  ... and {len(tokens) - 5} more")
    
    # Test satisfaction check
    print("\n=== Satisfaction Check ===")
    detected = {"ArBr_present": True, "ArB(OH)2_present": True}
    satisfied, matched = check_features_satisfy_reactants(detected, ["ArX*", "ArB*"])
    print(f"Suzuki reactants satisfied: {satisfied}")
    print(f"Matched: {matched}")
