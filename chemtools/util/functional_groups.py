"""
Functional Group Detection sourced from calculable_features.json.

All SMARTS definitions and metadata now live in the central
chemtools/featurizers/calculable_features.json specification, ensuring a
single source of truth for functional group logic.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .rdkit_helpers import rdkit_available, parse_smiles
from .smarts_cache import compile_smarts


@dataclass(frozen=True)
class FunctionalGroupDef:
    """Structured functional group information sourced from the feature spec."""

    name: str
    label: str
    smarts: Tuple[str, ...]
    text_patterns: Tuple[str, ...]
    category_tags: Tuple[str, ...]


_CATEGORY_PREFERRED_ORDER = [
    "oxygen",
    "nitrogen",
    "sulfur",
    "phosphorus",
    "halides",
    "aromatic",
    "unsaturated",
    "protecting_groups",
    "leaving_groups",
]
_UNCATEGORIZED_TAG = "other"
_SPEC_PATH = Path(__file__).resolve().parents[1] / "featurizers" / "calculable_features.json"


@lru_cache(maxsize=1)
def _load_feature_spec() -> Dict[str, object]:
    """Load the shared calculable feature specification from disk."""
    with _SPEC_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@lru_cache(maxsize=1)
def _load_group_definitions() -> Dict[str, FunctionalGroupDef]:
    """Load functional group definitions from the shared feature specification."""
    spec = _load_feature_spec()
    entries = spec.get("functional_groups") or []
    defs: Dict[str, FunctionalGroupDef] = {}

    for entry in entries:
        name = entry.get("name")
        detect = entry.get("detect") or {}
        smarts = tuple(detect.get("smarts_any") or [])
        if not name or not smarts:
            # Skip malformed entries rather than crashing downstream code.
            continue
        label = entry.get("label") or name.replace("_", " ").title()
        text_patterns = tuple(entry.get("text_patterns") or [])
        category_tags = tuple(entry.get("category_tags") or [])
        defs[name] = FunctionalGroupDef(
            name=name,
            label=label,
            smarts=smarts,
            text_patterns=text_patterns,
            category_tags=category_tags,
        )

    return defs


def get_group_definition(name: str) -> Optional[FunctionalGroupDef]:
    """Return the FunctionalGroupDef for a given canonical name."""
    return _load_group_definitions().get(name)


def iter_group_definitions() -> Iterable[FunctionalGroupDef]:
    """Yield all functional group definitions."""
    return _load_group_definitions().values()


@lru_cache(maxsize=1)
def _ordered_category_keys() -> List[str]:
    """Derive the ordered list of category tags based on spec metadata."""
    defs = _load_group_definitions()
    tags = set()
    for definition in defs.values():
        tags.update(definition.category_tags)

    ordered: List[str] = [tag for tag in _CATEGORY_PREFERRED_ORDER if tag in tags]
    extras = sorted(tags - set(_CATEGORY_PREFERRED_ORDER))
    ordered.extend(extras)
    ordered.append(_UNCATEGORIZED_TAG)
    return ordered


def _default_result_map() -> Dict[str, bool]:
    """Return a default result dict with all functional groups set to False."""
    return {name: False for name in _load_group_definitions().keys()}


def _has_smarts_match(mol, smarts_patterns: Tuple[str, ...]) -> bool:
    """Check whether any SMARTS pattern matches the RDKit molecule."""
    for smarts in smarts_patterns:
        pattern = compile_smarts(smarts, validate=False)
        if pattern is None:
            continue
        try:
            if mol.HasSubstructMatch(pattern):
                return True
        except Exception:
            continue
    return False


def _count_smarts_matches(mol, smarts_patterns: Tuple[str, ...]) -> int:
    """Count total matches for the provided SMARTS patterns."""
    total = 0
    for smarts in smarts_patterns:
        pattern = compile_smarts(smarts, validate=False)
        if pattern is None:
            continue
        try:
            matches = mol.GetSubstructMatches(pattern)
            total += len(matches)
        except Exception:
            continue
    return total


def _detect_with_text(smiles: str) -> Dict[str, bool]:
    """Fallback detection using simple substring checks."""
    lower = (smiles or "").lower()
    defs = _load_group_definitions()
    results: Dict[str, bool] = {}
    for name, definition in defs.items():
        if not definition.text_patterns:
            results[name] = False
            continue
        results[name] = any(pattern in lower for pattern in definition.text_patterns)
    return results


def detect_all(smiles: Optional[str]) -> Dict[str, bool]:
    """
    Detect all functional groups for a SMILES string.

    Args:
        smiles: Molecule SMILES.

    Returns:
        Mapping of functional group name to detection boolean.
    """
    if not smiles:
        return _default_result_map()

    if not rdkit_available():
        return _detect_with_text(smiles)

    mol = parse_smiles(smiles)
    if mol is None:
        return _detect_with_text(smiles)

    defs = _load_group_definitions()
    results: Dict[str, bool] = {}
    for name, definition in defs.items():
        results[name] = _has_smarts_match(mol, definition.smarts)
    return results


def get_functional_groups(smiles: Optional[str]) -> List[str]:
    """Return the list of functional groups present in the molecule."""
    detections = detect_all(smiles)
    return [name for name, present in detections.items() if present]


def has_functional_group(smiles: Optional[str], group_name: str) -> bool:
    """Check whether the molecule contains the specified functional group."""
    defs = _load_group_definitions()
    if group_name not in defs:
        return False
    result = detect_all(smiles)
    return result.get(group_name, False)


def count_functional_groups(smiles: Optional[str], group_name: str) -> int:
    """Count occurrences of the specified functional group."""
    defs = _load_group_definitions()
    definition = defs.get(group_name)
    if not smiles or definition is None or not rdkit_available():
        return 0
    mol = parse_smiles(smiles)
    if mol is None:
        return 0
    return _count_smarts_matches(mol, definition.smarts)


def get_group_categories(smiles: Optional[str]) -> Dict[str, List[str]]:
    """
    Categorize detected functional groups based on metadata tags.

    Returns:
        Dictionary mapping category tag -> list of group names within that category.
    """
    defs = _load_group_definitions()
    detections = get_functional_groups(smiles)
    categories = {tag: [] for tag in _ordered_category_keys()}

    for group in detections:
        definition = defs.get(group)
        if not definition:
            continue
        tags = definition.category_tags or ()
        if not tags:
            categories[_UNCATEGORIZED_TAG].append(group)
            continue
        for tag in tags:
            categories.setdefault(tag, []).append(group)

    return categories


def summarize_functional_groups(smiles: Optional[str]) -> str:
    """Human-readable summary of detected functional group categories."""
    categories = get_group_categories(smiles)
    lines: List[str] = []
    for tag in _ordered_category_keys():
        groups = categories.get(tag) or []
        if not groups:
            continue
        label = tag.replace("_", " ").title()
        lines.append(f"{label}: {', '.join(sorted(groups))}")
    return "\n".join(lines) if lines else "No functional groups detected"


def has_free_alcohol(smiles: Optional[str]) -> bool:
    """Check for a free (non-acid) alcohol."""
    if not smiles:
        return False
    has_carboxylic = has_functional_group(smiles, "carboxylic_acid")
    has_phenol_group = has_functional_group(smiles, "phenol")
    has_alcohol_group = has_functional_group(smiles, "alcohol")

    if has_carboxylic:
        count_oh = count_functional_groups(smiles, "alcohol")
        count_cooh = count_functional_groups(smiles, "carboxylic_acid")
        return count_oh > count_cooh

    return has_alcohol_group or has_phenol_group


def has_phenol(smiles: Optional[str]) -> bool:
    """Return True if a phenol functional group is present."""
    return has_functional_group(smiles, "phenol")


def has_sulfonamide(smiles: Optional[str]) -> bool:
    """Return True if a sulfonamide functional group is present."""
    return has_functional_group(smiles, "sulfonamide")


def has_hydroxylamine(smiles: Optional[str]) -> bool:
    """Return True if a hydroxylamine functional group is present."""
    return has_functional_group(smiles, "hydroxylamine")
