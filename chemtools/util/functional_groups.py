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
from types import MappingProxyType
from typing import Dict, Iterable, List, Optional, Tuple, Mapping

from .rdkit_helpers import rdkit_available, parse_smiles
from .smarts_cache import compile_smarts


@dataclass(frozen=True)
class FunctionalGroupDef:
    """Structured functional group information sourced from the feature spec."""

    name: str
    token: str
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
_DETECTION_CACHE_SIZE = 4096

# Fallback SMARTS for legacy group names still referenced by router/tests
_LEGACY_GROUP_FALLBACKS: Dict[str, Tuple[str, ...]] = {
    "acid": ("C(=O)[OH]",),
    "acyl_halide": ("[CX3](=O)[Cl,Br]",),
    "alcohol": ("[CX4][OX2H]",),
    "aldehyde": ("[CX3H](=O)",),
    "alkene": ("C=C",),
    "alkoxide": ("[O-][C,H]",),
    "alkyl_halide": ("[CX4][Cl,Br,I]",),
    "alpha_beta_unsaturated": ("C=C-[CX3]=O",),
    "aryl_halide": ("[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]",),
    "borane": ("[BH3,BH2,BH]",),
    "boron": ("[BX3;$(B(O)O),$(B(O)O),$(B(O)O)]",),
    "carbonyl": ("[CX3]=O",),
    "conjugated_diene": ("C=C-C=C",),
    "diene": ("C=C",),  # General diene pattern - molecule with two or more alkene groups
    "cyanide": ("[C-]#N",),
    "ester": ("[CX3](=O)[OX2][C,H]",),
    "grignard": ("[C,c][Mg][Br,Cl,I]",),
    "iodide": ("[I-]",),
    "ketone": ("[CX3](=O)[C]",),
    "nucleophile_n": ("[NX3;H1,H2]", "[nH]"),
    "nucleophile_o": ("[OX2H]",),
    "nucleophile_s": ("[SX2H]",),
    "organolithium": ("[C,c][Li]",),
    "organozinc": ("[C,c][Zn][Br,Cl,I]",),
    "phenol": ("[OH;$([OH]c1ccccc1)]",),  # Phenol pattern
    "terminal_alkene": ("C=C",),  # Terminal alkene pattern
    "terminal_alkyne": ("C#C[H]", "[C;H]#C"),
    "triflate": ("OS(=O)(=O)C(F)(F)F",),
    "vinyl_halide": ("C=C[Cl,Br,I]",),
}


@lru_cache(maxsize=1)
def _load_feature_spec() -> Dict[str, object]:
    """Load the shared calculable feature specification from disk."""
    with _SPEC_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@lru_cache(maxsize=1)
def _load_group_definitions() -> Dict[str, FunctionalGroupDef]:
    """Load functional group definitions from the shared feature specification."""
    spec = _load_feature_spec()
    defs: Dict[str, FunctionalGroupDef] = {}

    def register(entry: Dict[str, object]) -> None:
        detect = entry.get("detect") or {}
        smarts_list = detect.get("smarts_any") or []
        if isinstance(smarts_list, str):
            smarts_list = [smarts_list]
        smarts = tuple(smarts_list or [])
        if not smarts:
            return

        token = entry.get("token")
        name = entry.get("name")
        if not name:
            if isinstance(token, str) and token.endswith("_present"):
                name = token[:-8]
            elif isinstance(token, str):
                name = token
        if not isinstance(name, str) or not name:
            return

        if name in defs:
            return  # Prefer explicitly-declared functional_groups entries

        label = entry.get("label") or entry.get("why") or name.replace("_", " ").title()
        text_patterns = tuple(entry.get("text_patterns") or detect.get("text_patterns") or [])
        category_tags = entry.get("category_tags")
        if not category_tags:
            cat = entry.get("category")
            if cat:
                category_tags = [cat]
        category_tags = tuple(category_tags or [])
        token_value = token or f"{name}_present"

        defs[name] = FunctionalGroupDef(
            name=name,
            token=token_value,
            label=str(label),
            smarts=smarts,
            text_patterns=text_patterns,
            category_tags=category_tags,
        )

    for entry in spec.get("functional_groups") or []:
        register(entry)

    for feature in spec.get("features") or []:
        if feature.get("type") != "bool":
            continue
        if not isinstance(feature.get("detect"), dict):
            continue
        register(feature)

    for legacy_name, patterns in _LEGACY_GROUP_FALLBACKS.items():
        if legacy_name in defs:
            continue
        smarts = tuple(patterns)
        defs[legacy_name] = FunctionalGroupDef(
            name=legacy_name,
            token=f"{legacy_name}_present",
            label=legacy_name.replace("_", " ").title(),
            smarts=smarts,
            text_patterns=(),
            category_tags=(),
        )

    return defs


@lru_cache(maxsize=1)
def _definitions_by_token() -> Dict[str, FunctionalGroupDef]:
    defs = _load_group_definitions()
    return {definition.token: definition for definition in defs.values()}


def _build_smarts_index() -> Dict[str, Tuple[str, ...]]:
    """Materialize a name/token -> SMARTS tuple mapping for downstream consumers."""
    defs = _load_group_definitions()
    mapping: Dict[str, Tuple[str, ...]] = {}
    for name, definition in defs.items():
        mapping[name] = definition.smarts
        mapping[definition.token] = definition.smarts
    return mapping


def get_group_definition(name: str) -> Optional[FunctionalGroupDef]:
    """Return the FunctionalGroupDef for a given canonical name."""
    defs = _load_group_definitions()
    if name in defs:
        return defs[name]
    return _definitions_by_token().get(name)


def _resolve_group_token(name: str) -> Optional[str]:
    """Map a legacy group name or token to the canonical token string."""
    if not name:
        return None
    defs = _load_group_definitions()
    if name in defs:
        return defs[name].token
    if name in _definitions_by_token():
        return name
    return None


def iter_group_definitions() -> Iterable[FunctionalGroupDef]:
    """Yield all functional group definitions."""
    return _load_group_definitions().values()


# Expose immutable SMARTS lookup for consumers that need raw patterns
FUNCTIONAL_GROUP_SMARTS: Mapping[str, Tuple[str, ...]] = MappingProxyType(_build_smarts_index())


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
    return {definition.token: False for definition in _definitions_by_token().values()}


def _normalized_smiles(smiles: Optional[str]) -> str:
    """Normalize SMILES input for caching (empty string represents None)."""
    return (smiles or "").strip()


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
            results[definition.token] = False
            continue
        results[definition.token] = any(pattern in lower for pattern in definition.text_patterns)
    return results


def _detect_all_impl(smiles: str) -> Dict[str, bool]:
    """Internal detection implementation without caching."""
    if not smiles:
        return _default_result_map()

    if not rdkit_available():
        return _detect_with_text(smiles)

    mol = parse_smiles(smiles)
    if mol is None:
        return _detect_with_text(smiles)

    defs = _load_group_definitions()
    return {
        definition.token: _has_smarts_match(mol, definition.smarts)
        for definition in defs.values()
    }


@lru_cache(maxsize=_DETECTION_CACHE_SIZE)
def _detect_all_cached(smiles: str) -> Tuple[Tuple[str, bool], ...]:
    """LRU-cached detection results keyed by normalized SMILES."""
    return tuple(_detect_all_impl(smiles).items())


def detect_all(smiles: Optional[str]) -> Dict[str, bool]:
    """
    Detect all functional groups for a SMILES string.

    Args:
        smiles: Molecule SMILES.

    Returns:
        Mapping of functional-group tokens (``<name>_present``) to booleans.
    """
    key = _normalized_smiles(smiles)
    if not key:
        return _default_result_map()
    return dict(_detect_all_cached(key))


def detect_any(smiles_list: Iterable[Optional[str]], *, group_subset: Optional[Iterable[str]] = None) -> Dict[str, bool]:
    """
    Detect functional groups across multiple SMILES strings (logical OR of detections).

    Args:
        smiles_list: Iterable of SMILES strings.
        group_subset: Optional subset of group names to report (improves performance).

    Returns:
        Mapping of group name to boolean indicating presence in ANY provided SMILES.
    """
    if group_subset:
        ordered_subset: List[str] = []
        for name in group_subset:
            token = _resolve_group_token(name)
            if token and token not in ordered_subset:
                ordered_subset.append(token)
        aggregate = {name: False for name in ordered_subset}
    else:
        aggregate = _default_result_map()

    for smiles in smiles_list:
        key = _normalized_smiles(smiles)
        if not key:
            continue
        detections = detect_all(key)
        if group_subset:
            for name in ordered_subset:
                if aggregate.get(name):
                    continue
                if detections.get(name):
                    aggregate[name] = True
        else:
            for name, present in detections.items():
                if present:
                    aggregate[name] = True
    return aggregate


def get_functional_groups(smiles: Optional[str]) -> List[str]:
    """Return the list of functional groups present in the molecule."""
    detections = detect_all(smiles)
    return [name for name, present in detections.items() if present]


def has_functional_group(smiles: Optional[str], group_name: str) -> bool:
    """Check whether the molecule contains the specified functional group."""
    token = _resolve_group_token(group_name)
    if not token:
        return False
    result = detect_all(smiles)
    return result.get(token, False)


def count_functional_groups(smiles: Optional[str], group_name: str) -> int:
    """Count occurrences of the specified functional group."""
    defs = _load_group_definitions()
    definition = defs.get(group_name)
    if definition is None:
        definition = _definitions_by_token().get(group_name)
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
    defs_by_token = _definitions_by_token()
    detections = get_functional_groups(smiles)
    categories = {tag: [] for tag in _ordered_category_keys()}

    for group in detections:
        definition = defs_by_token.get(group)
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
