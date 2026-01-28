"""Reactant classification using the unified feature system.

This module provides SMARTS-driven reactant classification that delegates
to the unified feature detection system (chemtools.featurizers.calculable).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

from .taxonomy import _load_reactant_types_raw


def _get_calculable() -> Any:
    """Lazy import of the calculable module to avoid circular imports."""
    from ... import calculable as _calc
    return _calc


@dataclass(frozen=True)
class ReactantMatch:
    """Structured result for a SMARTS classification."""

    category: str
    member_type: str
    name: str
    group: str
    smarts: str
    category_smarts: Optional[str]
    description: str
    specificity: int
    is_general: bool


def reactant_category_for(match: ReactantMatch) -> str:
    """Return the canonical category associated with a match."""
    return match.category


def build_reactant_lookup() -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (alias -> canonical reactant id, canonical id -> category id) mappings."""
    definitions = _load_reactant_types_raw()

    alias_map: Dict[str, str] = {}
    id_to_category: Dict[str, str] = {}

    for category_id, data in definitions.items():
        for member in data.get("members", []):
            canonical = member.get("id")
            if not canonical:
                continue
            id_to_category[canonical] = category_id
            alias_map[canonical.lower()] = canonical
            for alias in member.get("aliases", []):
                if alias:
                    alias_map[alias.lower()] = canonical
    return alias_map, id_to_category


def iter_reactant_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[ReactantMatch]:
    """Return all SMARTS matches for ``smiles``.

    Now uses the unified feature system for detection, returning all matching
    reactant types as ReactantMatch objects.

    Args:
        smiles: SMILES string to analyze
        reactant_types: Deprecated, no longer used (for backward compatibility only)

    Returns:
        List of ReactantMatch objects for all detected reactant types
    """
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored. "
            "Classification now uses the unified feature system.",
            DeprecationWarning,
            stacklevel=2,
        )

    smiles = (smiles or "").strip()
    if not smiles:
        return []

    # Use the unified feature system to get all reactant type features
    _calc = _get_calculable()
    reactant_features = _calc.get_reactant_type_features(smiles)
    if not reactant_features:
        return []

    matches: List[ReactantMatch] = []

    # Get all detected member types
    member_types = reactant_features.get("member_types", [])

    # Load feature specs to get metadata
    spec = _calc._load_feature_spec()
    features_by_token = {f["token"]: f for f in spec.get("features", [])}

    # Build matches from detected member-level features
    for member_id in member_types:
        token = f"{member_id}_reactant"
        if not reactant_features.get(token, False):
            continue

        feature_spec = features_by_token.get(token, {})
        metadata = feature_spec.get("metadata", {})

        # Extract SMARTS pattern
        detect = feature_spec.get("detect", {})
        smarts_patterns = detect.get("smarts_any", [])
        smarts = smarts_patterns[0] if smarts_patterns else ""

        category = metadata.get("reactant_category", "")
        priority = metadata.get("priority", 1)

        matches.append(
            ReactantMatch(
                category=category,
                member_type=member_id,
                name=metadata.get("member_name", metadata.get("name", "")),
                group=metadata.get("group", ""),
                smarts=smarts,
                category_smarts=metadata.get("category_smarts"),
                description=metadata.get("description", ""),
                specificity=priority,  # Use priority as specificity
                is_general=priority < 2,  # Heuristic: priority 1 is general
            )
        )

    # Sort by (priority descending, member id for determinism)
    matches.sort(key=lambda m: (-m.specificity, m.member_type))

    return matches


def classify_reactant_smiles(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[ReactantMatch]:
    """Return the most specific SMARTS match for the SMILES input.

    Now delegates to the unified feature system (chemtools.featurizers.calculable)
    for reactant detection, but maintains backward compatibility by returning
    a ReactantMatch with the expected structure.

    Args:
        smiles: SMILES string to classify
        reactant_types: Deprecated, no longer used (for backward compatibility only)

    Returns:
        ReactantMatch with category, member_type, name, etc., or None if no match
    """
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored. "
            "Classification now uses the unified feature system.",
            DeprecationWarning,
            stacklevel=2,
        )

    # Use the unified feature system
    _calc = _get_calculable()
    result = _calc.classify_reactant_smiles(smiles)
    if result is None:
        return None

    # Convert to ReactantMatch structure expected by this module
    return ReactantMatch(
        category=result.get("category", ""),
        member_type=result.get("member_type", ""),
        name=result.get("name", ""),
        group=result.get("group", ""),
        smarts=result.get("smarts", ""),
        category_smarts=result.get("category_smarts"),
        description=result.get("description", ""),
        specificity=result.get("specificity", 0),
        is_general=result.get("is_general", False),
    )


def classify_reactant_category(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[str]:
    """Convenience shortcut returning only the category id."""
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    best = classify_reactant_smiles(smiles)
    return best.category if best else None


def classify_reactant_group(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[str]:
    """Convenience shortcut returning the functional group label."""
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    best = classify_reactant_smiles(smiles)
    return best.group if best else None


def classify_reactant_batch(
    smiles_list: Iterable[str], reactant_types: Optional[Dict[str, dict]] = None
) -> List[Optional[ReactantMatch]]:
    """Batch classification wrapper for reactant type detection."""
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    return [classify_reactant_smiles(smiles) for smiles in smiles_list]


def get_reactant_category_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[str]:
    """Return the set of categories matched by the SMARTS hierarchy."""
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    # Use unified feature system directly
    _calc = _get_calculable()
    reactant_features = _calc.get_reactant_type_features(smiles)
    if not reactant_features:
        return []
    categories = reactant_features.get("categories", [])
    if not categories:
        return []

    index = _calc._reactant_feature_index()
    priorities: Dict[str, int] = {}
    for member in index.get("members", []):
        token = member.get("token")
        category_id = member.get("category_id")
        if not token or not category_id:
            continue
        if reactant_features.get(token, False):
            priority = int(member.get("priority", 1))
            if priority > priorities.get(category_id, -1):
                priorities[category_id] = priority

    for category_id, entry in index.get("categories", {}).items():
        token = entry.get("token")
        if (
            token
            and reactant_features.get(token, False)
            and category_id not in priorities
        ):
            priorities[category_id] = 0

    if not priorities:
        return sorted(categories)

    return [
        cid for cid, _ in sorted(priorities.items(), key=lambda kv: (-kv[1], kv[0]))
    ]


def get_all_reactant_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[ReactantMatch]:
    """Alias retained for backwards compatibility with HTE scripts."""
    if reactant_types is not None:
        warnings.warn(
            "The 'reactant_types' parameter is deprecated and will be ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    return iter_reactant_matches(smiles)


__all__ = [
    "ReactantMatch",
    "reactant_category_for",
    "build_reactant_lookup",
    "iter_reactant_matches",
    "classify_reactant_smiles",
    "classify_reactant_category",
    "classify_reactant_group",
    "classify_reactant_batch",
    "get_reactant_category_matches",
    "get_all_reactant_matches",
]
