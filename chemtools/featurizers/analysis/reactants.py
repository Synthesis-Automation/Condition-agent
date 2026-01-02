"""
Reactant and reaction type utilities shared across ChemTools and HTE pipelines.

This module centralises the SMARTS-driven reactant classification helpers using
the unified feature detection system (``chemtools.featurizers.calculable``) as
the single source of truth for reactant identification (organic_compounds).

Note: This module now delegates reactant classification to the unified feature
system while maintaining backward compatibility with the legacy API.
"""

from __future__ import annotations

import copy
import json
import re
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .. import calculable as _calculable

from ...taxonomy import reaction_catalog as _reaction_catalog
from ...util import rdkit_helpers
from ...util.smarts_cache import compile_smarts
from ._registry import clear_registry_cache, get_registry

TAXONOMY_DATA_DIR = Path(__file__).resolve().parents[2] / "taxonomy" / "data"
REACTANT_TYPES_FILE = TAXONOMY_DATA_DIR / "organic_compounds.v1.3.json"
REACTION_TYPES_FILE = _reaction_catalog.REACTION_TYPES_FILE

# General categories that should be deprioritised when picking the "best" match.
GENERAL_REACTANT_CATEGORIES = {"R-H", "Ar-H"}


def _get_calculable() -> Any:
    """Lazy import of the calculable module to avoid circular imports."""
    from .. import calculable as _calc
    return _calc

# Manual overrides sourced from the original z-score pipeline.
CSV_REACTANT_OVERRIDES: Dict[str, str] = {
    "ArBr": "Ar-Br",
    "ArCl": "Ar-Cl",
    "ArI": "Ar-I",
    "ArF": "Ar-F",
    "ArOH": "Ar-OH",
    "ArOSO2R": "Ar-OTs",
    "Alkyl-Br": "R-Br",
    "Alkyl-Cl": "R-Cl",
    "Alkyl-I": "R-I",
    "Alkyl-OSO2R": "R-OSO2R",
    "alkene-Br": "Vinyl-Br",
    "alkene-I": "Vinyl-I",
    "ArB(OR)2": "Ar-B(OR)2",
    "ArB(OH)2": "Ar-B(OH)2",
    "ArBF3K": "Ar-BF3K",
    "RNH2 a-branch": "Any-NH2",
    "RNH2": "Any-NH2",
    "R2NH a-branch": "Any-NHR",
    "R2NH": "Any-NHR",
    "ArNH2": "Ar-NH2",
    "ArNHR": "Ar-NHR",
    "Ar2NH": "Ar-NHAr",
    "arom. NH": "AromN-H",
    "alkeneB(OR)2": "Vinyl-B(OR)2",
    "Alkyl-B(OH)2": "R-B(OH)2",
    "Alkyl-B(OR)2": "Any-B(OR)2",
    "Alkyl-BF3K": "R-BF3K",
    "Alkyl-OH a-branch": "R2CH-OH",
    "Alkyl-OH primary": "RCH2-OH",
    "Lactam": "Any-Lactam",
    "RCONH2": "Any-CONHR",
    "ROC(O)NR2": "Any-Carbamate",
    "Urea": "Any-Urea",
    "Alkyl-M": "R-M",
    "Alkyl-H acidic": "R-H-acidic",
    "Alkyl-H": "R-H",
    "alkene": "Any-Alkene",
    "alkyne": "Any-Alkyne",
    "Ar-H": "Ar-H",
    "enolether": "Any-EnolEther",
    "RSH": "Any-SH",
    "RCO2H or M": "Any-CO2H",
}

CSV_REACTION_OVERRIDES: Dict[str, str] = {
    "Buchwald-Hartwig": "Buchwald-Hartwig",
    "Suzuki-Miyaura": "Suzuki-Miyaura",
    "Suzuki-Miyaura, in situ": "Suzuki-Miyaura-in-situ",
    "Arylation, acidic C-H": "Arylation-acidic-C-H",
    "Amide coupling": "Amide-coupling",
    "CN-Coupling": "CN-Coupling",
    "CO-Coupling": "CO-Coupling",
    "Condensation": "Condensation",
    "CH-Activation": "CH-Activation",
    "Negishi, in-situ": "Negishi-in-situ",
    "Cyclization": "Cyclization",
    "Borylation, Miyaura": "Borylation-Miyaura",
    "Alkylation": "Alkylation",
    "Deprotection": "Deprotection",
    "Negishi": "Negishi",
    "Heck": "Heck",
    "CC-Coupling": "CC-Coupling",
    "SNAr": "SNAr",
    "Hydrolysis": "Hydrolysis",
    "Salt formation": "Salt-formation",
    "Sonogashira": "Sonogashira",
    "Stetter": "Stetter",
    "Cyanation": "Cyanation",
    "Oxidation": "Oxidation",
    "Activation": "Activation",
    "Hydrodehalogenation": "Hydrodehalogenation",
    "Glycosidation": "Glycosidation",
    "Stannylation": "Stannylation",
    "Dehydration": "Dehydration",
    "Dimerization, reductive": "Dimerization-reductive",
    "Mitsunobu": "Mitsunobu",
    "CS-Coupling": "CS-Coupling",
    "Borylation, C-H": "Borylation-C-H",
    "Wittig": "Wittig",
    "Sandmeyer": "Sandmeyer",
    "Addition": "Addition",
    "Hydration": "Hydration",
    "Reduction": "Reduction",
    "Deoxyfluorination": "Deoxyfluorination",
    "Chlorination": "Halogenation-Chlorination",
    "Fluorination, oxidative": "Fluorination-oxidative",
    "Protection": "Protection",
}


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


def _pick_legacy_alias(entity_type: str, entity_id: str) -> Optional[str]:
    registry = get_registry()
    if registry is None:
        return None

    matches = [
        record
        for record in registry.aliases.values()
        if record.entity_type == entity_type and record.entity_id == entity_id
    ]
    if not matches:
        return None

    if entity_type == "reactant_type":
        filtered = [
            record
            for record in matches
            if (record.notes or "").lower() not in {"member_id", "member"}
        ]
        if filtered:
            matches = filtered

    def _is_preferred(alias: str) -> bool:
        if not alias:
            return False
        if alias.lower() != alias:
            return True
        return any(ch in alias for ch in ("*", "-", " "))

    # Prefer aliases explicitly marked as original/legacy.
    for record in matches:
        note = (record.notes or "").lower()
        if note in {"original_id", "legacy_id"} and record.alias:
            if record.alias.lower() != entity_id.lower():
                return record.alias

    for record in matches:
        if record.alias and record.alias.lower() != entity_id.lower() and _is_preferred(record.alias):
            return record.alias

    for record in matches:
        if record.alias:
            return record.alias

    return None


@lru_cache(maxsize=1)
def _load_reactant_types_raw() -> Dict[str, dict]:
    registry = get_registry()
    if registry is None:
        return _load_reactant_types_from_file(REACTANT_TYPES_FILE)

    definitions: Dict[str, dict] = {}
    for reactant_id, reactant in registry.reactant_types.items():
        legacy_id = _pick_legacy_alias("reactant_type", reactant_id)
        members: List[dict] = []
        for member in reactant.members:
            members.append(
                {
                    "id": member.id,
                    "name": member.name,
                    "smarts": member.smarts,
                    "aliases": list(member.aliases),
                    "metadata": copy.deepcopy(member.metadata),
                }
            )
        definitions[reactant_id] = {
            "id": reactant_id,
            "name": reactant.name,
            "description": reactant.description,
            "category": reactant.category,
            "smarts": reactant.smarts,
            "group": (reactant.metadata or {}).get("group", ""),
            "aliases": list(reactant.aliases),
            "metadata": copy.deepcopy(reactant.metadata),
            "members": members,
            "legacy_id": legacy_id,
        }
    return definitions


def _load_reactant_types_from_file(path: Path) -> Dict[str, dict]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if isinstance(payload, dict) and "compounds" in payload:
        return _load_reactant_types_from_compounds(payload.get("compounds") or [])

    if isinstance(payload, dict):
        entries = payload.get("entries") or []
    else:
        entries = payload

    definitions: Dict[str, dict] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        entry_id = entry.get("id")
        if not entry_id:
            continue
        entry_meta = dict(entry.get("metadata") or {})
        feature_token = entry.get("feature_token")
        if feature_token and "feature_token" not in entry_meta:
            entry_meta["feature_token"] = feature_token

        members: List[dict] = []
        for member in entry.get("members", []):
            if not isinstance(member, dict):
                continue
            member_id = member.get("id")
            if not member_id:
                continue
            member_meta = dict(member.get("metadata") or {})
            member_token = member.get("feature_token")
            if member_token and "feature_token" not in member_meta:
                member_meta["feature_token"] = member_token
            members.append(
                {
                    "id": member_id,
                    "name": member.get("name", member_id),
                    "smarts": member.get("smarts"),
                    "aliases": list(member.get("aliases", [])),
                    "metadata": member_meta,
                }
            )

        definitions[entry_id] = {
            "id": entry_id,
            "name": entry.get("name", entry_id),
            "description": entry.get("description"),
            "category": entry.get("category"),
            "smarts": entry.get("smarts"),
            "group": (entry_meta or {}).get("group", ""),
            "aliases": list(entry.get("aliases", [])),
            "metadata": entry_meta,
            "members": members,
            "legacy_id": entry.get("legacy_id") or entry_id,
        }

    return definitions


def _load_reactant_types_from_compounds(compounds: Iterable[dict]) -> Dict[str, dict]:
    definitions: Dict[str, dict] = {}
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        compound_id = entry.get("id")
        if not compound_id:
            continue
        name = entry.get("name", compound_id)
        aliases = list(entry.get("aliases") or [])
        smarts_any = entry.get("smarts_any") or entry.get("smarts")
        smarts_list: List[str] = []
        if isinstance(smarts_any, str):
            smarts_list = [smarts_any]
        elif isinstance(smarts_any, list):
            smarts_list = [s for s in smarts_any if isinstance(s, str)]
        smarts = smarts_list[0] if smarts_list else ""
        group = str(entry.get("B") or "")

        definitions[compound_id] = {
            "id": compound_id,
            "name": name,
            "description": entry.get("description"),
            "category": compound_id,
            "smarts": smarts,
            "group": group,
            "aliases": aliases,
            "metadata": {"group": group},
            "members": [
                {
                    "id": compound_id,
                    "name": name,
                    "smarts": smarts,
                    "aliases": aliases,
                    "metadata": {"legacy_taxonomy_id": compound_id, "group": group},
                }
            ],
            "legacy_id": compound_id,
        }

    return definitions


def get_reactant_type_definitions() -> Dict[str, dict]:
    """Return a deep copy of the reactant type taxonomy (organic_compounds-backed)."""
    return copy.deepcopy(_load_reactant_types_raw())


def get_reactant_types_file() -> Path:
    """Return the path to the canonical reactant type JSON file (organic_compounds)."""
    return REACTANT_TYPES_FILE


def clear_reactant_type_cache() -> None:
    """Clear the cached reactant type definitions (useful after editing JSON)."""
    _load_reactant_types_raw.cache_clear()
    _reactant_alias_index.cache_clear()
    clear_registry_cache()


@lru_cache(maxsize=None)
def _reactant_alias_index() -> Dict[str, str]:
    """Build a lookup of aliases → canonical category IDs."""
    definitions = _load_reactant_types_raw()
    alias_map: Dict[str, str] = {}
    for category, data in definitions.items():
        legacy_category = data.get("legacy_id") or category
        alias_map[category.lower()] = legacy_category
        alias_map[legacy_category.lower()] = legacy_category
        for alias in data.get("aliases", []):
            alias_map[alias.lower()] = legacy_category
        for member in data.get("members", []):
            member_id = member.get("id", "")
            if member_id:
                alias_map[member_id.lower()] = legacy_category
            for alias in member.get("aliases", []):
                alias_map[alias.lower()] = legacy_category
    alias_map.update({k.lower(): v for k, v in CSV_REACTANT_OVERRIDES.items()})
    return alias_map


def normalize_reactant_identifier(label: str) -> Optional[str]:
    """Return the canonical reactant category id for ``label``."""
    if not label:
        return None
    return _reactant_alias_index().get(label.strip().lower())


def build_reactant_lookup() -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Return (alias -> canonical reactant id, canonical id -> category id) mappings.
    """
    definitions = _load_reactant_types_raw()

    alias_map: Dict[str, str] = {}
    id_to_category: Dict[str, str] = {}

    for category_id, data in definitions.items():
        legacy_category = data.get("legacy_id") or category_id
        for member in data.get("members", []):
            canonical = member.get("id")
            if not canonical:
                continue
            id_to_category[canonical] = legacy_category
            alias_map[canonical.lower()] = canonical
            for alias in member.get("aliases", []):
                if alias:
                    alias_map[alias.lower()] = canonical
    return alias_map, id_to_category


def reactant_category_for(match: ReactantMatch) -> str:
    """Return the canonical category associated with a match."""
    return match.category


def iter_reactant_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[ReactantMatch]:
    """
    Return all SMARTS matches for ``smiles``.
    
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
    categories = reactant_features.get("categories", [])
    
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
        is_general = category in GENERAL_REACTANT_CATEGORIES
        
        matches.append(
            ReactantMatch(
                category=category,
                member_type=member_id,
                name=metadata.get("member_name", metadata.get("name", "")),
                group=metadata.get("group", ""),
                smarts=smarts,
                category_smarts=metadata.get("category_smarts"),
                description=metadata.get("description", ""),
                specificity=len(smarts) if smarts else 0,
                is_general=is_general,
            )
        )
    
    # Sort by (general? -> False first, specificity descending, member id for determinism)
    matches.sort(key=lambda m: (m.is_general, -m.specificity, m.member_type))
    
    return matches


def classify_reactant_smiles(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[ReactantMatch]:
    """
    Return the most specific SMARTS match for the SMILES input.
    
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
    """Batch classification wrapper mirroring the legacy helper."""
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
    return sorted(reactant_features.get("categories", []))


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


@lru_cache(maxsize=1)
def _load_reaction_types_raw() -> Dict[str, dict]:
    definitions, _ = _reaction_catalog.load_reaction_catalog()
    categories: Dict[str, dict] = {}

    for reaction in definitions.values():
        category_id = reaction.category or "uncategorised"
        bucket = categories.setdefault(
            category_id,
            {
                "category_id": category_id,
                "category": category_id,
                "description": None,
                "reactions": [],
            },
        )
        bucket["reactions"].append(
            {
                "id": reaction.id,
                "canonical_id": reaction.id,
                "name": reaction.name,
                "aliases": list(reaction.aliases),
                "description": reaction.description,
                "reactants": copy.deepcopy(reaction.reactants),
                "metadata": copy.deepcopy(reaction.metadata),
                "catalysts": list(reaction.catalysts),
                "conditions": reaction.conditions,
                "reference_reactions": list(reaction.reference_reactions),
                "notes": reaction.notes,
            }
        )

    return categories


def get_reaction_type_definitions() -> Dict[str, dict]:
    """Return a deep copy of the reaction type taxonomy."""
    return copy.deepcopy(_load_reaction_types_raw())


def get_reaction_types_file() -> Path:
    """Return the path to the canonical reaction type JSON file."""
    return REACTION_TYPES_FILE


def clear_reaction_type_cache() -> None:
    """Clear the cached reaction type definitions."""
    _load_reaction_types_raw.cache_clear()
    _reaction_indices.cache_clear()
    clear_registry_cache()


@lru_cache(maxsize=1)
def _reaction_indices() -> Tuple[Dict[str, dict], Dict[str, str]]:
    """Build lookup tables (id -> metadata, alias -> id)."""
    definitions = _load_reaction_types_raw()
    id_to_meta: Dict[str, dict] = {}
    alias_to_id: Dict[str, str] = {}

    def register_alias(alias: str, reaction_id: str) -> None:
        if not alias:
            return
        alias_to_id[alias.lower()] = reaction_id

    for category_key, category_data in definitions.items():
        category_label = category_data.get("category")
        for reaction in category_data.get("reactions", []):
            reaction_id = reaction.get("id")
            if not reaction_id:
                continue

            metadata = copy.deepcopy(reaction)
            metadata["category_key"] = category_key
            metadata["category"] = category_label
            id_to_meta[reaction_id] = metadata

            register_alias(reaction_id, reaction_id)
            register_alias(reaction.get("name", ""), reaction_id)
            for alias in reaction.get("aliases", []):
                register_alias(alias, reaction_id)

    return id_to_meta, alias_to_id


def list_reaction_type_ids() -> List[str]:
    """Return all reaction ids defined in the taxonomy."""
    id_to_meta, _ = _reaction_indices()
    return sorted(id_to_meta.keys())


def describe_reaction_type(reaction_id: str) -> Optional[dict]:
    """Return metadata for the given reaction id (category, aliases, description, etc.)."""
    id_to_meta, _ = _reaction_indices()
    return copy.deepcopy(id_to_meta.get(reaction_id))


def normalize_reaction_type(label: str) -> Optional[str]:
    """Normalise a raw reaction descriptor (name/alias/id) to the canonical id."""
    if not label:
        return None
    label = label.strip()
    if not label:
        return None
    _, alias_to_id = _reaction_indices()
    key = label.lower()
    resolved = alias_to_id.get(key)
    if resolved:
        return resolved
    slug = re.sub(r"[^0-9a-z]+", "_", key).strip("_")
    if not slug:
        return None
    return alias_to_id.get(slug)


def build_reaction_lookup() -> Tuple[Dict[str, dict], Dict[str, str]]:
    """Expose the cached lookup tables (id -> meta, alias -> id)."""
    id_to_meta, alias_to_id = _reaction_indices()
    return copy.deepcopy(id_to_meta), copy.deepcopy(alias_to_id)


def iter_reactions_for_category(category_key: str) -> List[dict]:
    """Return all reactions defined under a given category key."""
    definitions = _load_reaction_types_raw()
    category = definitions.get(category_key)
    if not category:
        return []
    return copy.deepcopy(category.get("reactions", []))


def required_reactant_categories(reaction_id: str) -> Optional[Dict[str, List[str]]]:
    """
    Return the reactant slot requirements for a reaction id.

    The taxonomy v2 schema stores reactants as slot lists, e.g.
    ``{"electrophiles": [...], "nucleophiles": [...], "acids": [...]}``.
    """
    metadata = describe_reaction_type(reaction_id)
    if not metadata:
        return None
    reactants = metadata.get("reactants")
    if isinstance(reactants, dict):
        return copy.deepcopy(reactants)
    return None


__all__ = [
    "ReactantMatch",
    "GENERAL_REACTANT_CATEGORIES",
    "CSV_REACTANT_OVERRIDES",
    "CSV_REACTION_OVERRIDES",
    "get_reactant_type_definitions",
    "get_reactant_types_file",
    "clear_reactant_type_cache",
    "normalize_reactant_identifier",
    "build_reactant_lookup",
    "reactant_category_for",
    "iter_reactant_matches",
    "classify_reactant_smiles",
    "classify_reactant_category",
    "classify_reactant_group",
    "classify_reactant_batch",
    "get_reactant_category_matches",
    "get_all_reactant_matches",
    "get_reaction_type_definitions",
    "get_reaction_types_file",
    "clear_reaction_type_cache",
    "list_reaction_type_ids",
    "describe_reaction_type",
    "normalize_reaction_type",
    "build_reaction_lookup",
    "iter_reactions_for_category",
    "required_reactant_categories",
]
