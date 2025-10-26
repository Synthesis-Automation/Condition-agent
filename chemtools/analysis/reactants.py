"""
Reactant and reaction type utilities shared across ChemTools and HTE pipelines.

This module centralises the SMARTS-driven reactant classification helpers using
the unified taxonomy registry (``chemtools.taxonomy``) as the single source of
truth for reactant and reaction definitions.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from ..util import rdkit_helpers
from ._registry import clear_registry_cache, get_registry

TAXONOMY_DATA_DIR = Path(__file__).resolve().parent.parent / "taxonomy" / "data"
REACTANT_TYPES_FILE = TAXONOMY_DATA_DIR / "reactant_types.json"
REACTION_TYPES_FILE = TAXONOMY_DATA_DIR / "reaction_types.json"

# General categories that should be deprioritised when picking the "best" match.
GENERAL_REACTANT_CATEGORIES = {"Alkyl-C-H", "ArH"}

# Manual overrides sourced from the original z-score pipeline.
CSV_REACTANT_OVERRIDES: Dict[str, str] = {
    "ArBr": "ArBr",
    "ArCl": "ArCl",
    "ArI": "ArI",
    "ArF": "ArF",
    "ArOH": "ArOH",
    "ArOSO2R": "ArOSO2R",
    "Alkyl-Br": "Alkyl-Br",
    "Alkyl-Cl": "Alkyl-Cl",
    "Alkyl-I": "Alkyl-I",
    "Alkyl-OSO2R": "Alkyl-OSO2R",
    "alkene-Br": "alkene-Br",
    "alkene-I": "alkene-I",
    "ArB(OR)2": "ArB(OR)2",
    "ArB(OH)2": "ArB(OH)2",
    "ArBF3K": "ArBF3K",
    "RNH2 a-branch": "RNH2-a-branch",
    "RNH2": "RNH2",
    "R2NH a-branch": "R2NH-a-branch",
    "R2NH": "R2NH",
    "ArNH2": "ArNH2",
    "ArNHR": "ArNHR",
    "Ar2NH": "Ar2NH",
    "arom. NH": "arom-NH",
    "alkeneB(OR)2": "alkeneB(OR)2",
    "Alkyl-B(OH)2": "Alkyl-B(OH)2",
    "Alkyl-B(OR)2": "Alkyl-B(OR)2",
    "Alkyl-BF3K": "Alkyl-BF3K",
    "Alkyl-OH a-branch": "ROH-a-branch",
    "Alkyl-OH primary": "ROH-primary",
    "Lactam": "Lactam",
    "RCONH2": "RCONH2",
    "ROC(O)NR2": "Carbamate",
    "Urea": "Urea",
    "Alkyl-M": "Alkyl-M",
    "Alkyl-H acidic": "Alkyl-H-acidic",
    "Alkyl-H": "Alkyl-H",
    "alkene": "alkene",
    "alkyne": "alkyne",
    "Ar-H": "Ar-H",
    "enolether": "enol-ether",
    "RSH": "RSH",
    "RCO2H or M": "RCO2H",
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
        return {}

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


def get_reactant_type_definitions() -> Dict[str, dict]:
    """Return a deep copy of the reactant type taxonomy."""
    return copy.deepcopy(_load_reactant_types_raw())


def get_reactant_types_file() -> Path:
    """Return the path to the canonical reactant type JSON file."""
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
    """Return all SMARTS matches for ``smiles``."""
    smiles = (smiles or "").strip()
    if not smiles:
        return []

    mol = rdkit_helpers.parse_smiles(smiles)
    if mol is None:
        return []

    definitions = reactant_types if reactant_types is not None else _load_reactant_types_raw()
    matches: List[ReactantMatch] = []

    for category_id, data in definitions.items():
        legacy_category = data.get("legacy_id") or category_id
        category_smarts = data.get("smarts")
        category_description = data.get("description", "")
        group = data.get("group", "")
        is_general = legacy_category in GENERAL_REACTANT_CATEGORIES

        for member in data.get("members", []):
            smarts = member.get("smarts")
            if not smarts:
                continue

            pattern = rdkit_helpers.parse_smarts(smarts)
            if pattern is None:
                continue

            try:
                has_match = mol.HasSubstructMatch(pattern)
            except Exception:
                has_match = False

            if has_match:
                matches.append(
                    ReactantMatch(
                        category=legacy_category,
                        member_type=member.get("id", ""),
                        name=member.get("name", ""),
                        group=group,
                        smarts=smarts,
                        category_smarts=category_smarts,
                        description=category_description,
                        specificity=len(smarts),
                        is_general=is_general,
                    )
                )

    return matches


def classify_reactant_smiles(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[ReactantMatch]:
    """Return the most specific SMARTS match for the SMILES input."""
    matches = iter_reactant_matches(smiles, reactant_types)
    if not matches:
        return None

    # Sort by (general? -> False first, specificity descending, member id for determinism).
    matches.sort(key=lambda m: (m.is_general, -m.specificity, m.member_type))
    return matches[0]


def classify_reactant_category(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[str]:
    """Convenience shortcut returning only the category id."""
    best = classify_reactant_smiles(smiles, reactant_types)
    return best.category if best else None


def classify_reactant_group(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> Optional[str]:
    """Convenience shortcut returning the functional group label."""
    best = classify_reactant_smiles(smiles, reactant_types)
    return best.group if best else None


def classify_reactant_batch(
    smiles_list: Iterable[str], reactant_types: Optional[Dict[str, dict]] = None
) -> List[Optional[ReactantMatch]]:
    """Batch classification wrapper mirroring the legacy helper."""
    return [classify_reactant_smiles(smiles, reactant_types) for smiles in smiles_list]


def get_reactant_category_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[str]:
    """Return the set of categories matched by the SMARTS hierarchy."""
    matches = iter_reactant_matches(smiles, reactant_types)
    return sorted({match.category for match in matches})


def get_all_reactant_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[ReactantMatch]:
    """Alias retained for backwards compatibility with HTE scripts."""
    return iter_reactant_matches(smiles, reactant_types)


@lru_cache(maxsize=1)
def _load_reaction_types_raw() -> Dict[str, dict]:
    registry = get_registry()
    if registry is None:
        return {}

    categories: Dict[str, dict] = {}
    for category_id, category in registry.reaction_categories.items():
        categories[category_id] = {
            "category_id": category_id,
            "category": category.name or category_id,
            "description": category.description,
            "reactions": [],
        }

    for reaction in registry.reaction_types.values():
        category_id = reaction.category_id or "uncategorised"
        legacy_reaction_id = _pick_legacy_alias("reaction_type", reaction.id) or reaction.id
        bucket = categories.setdefault(
            category_id,
            {
                "category_id": category_id,
                "category": (
                    registry.reaction_categories.get(category_id).name
                    if registry and category_id in registry.reaction_categories
                    else category_id
                ),
                "description": (
                    registry.reaction_categories.get(category_id).description
                    if registry and category_id in registry.reaction_categories
                    else None
                ),
                "reactions": [],
            },
        )

        reactant_blocks: List[List[str]] = []
        for requirement in reaction.reactants:
            tokens = list(requirement.original_tokens or [])
            if not tokens and requirement.reactant_type_id:
                tokens = [requirement.reactant_type_id]

            normalized_tokens: List[str] = []
            for token in tokens:
                alias = _pick_legacy_alias("reactant_type", token)
                if alias:
                    normalized_tokens.append(alias)
                else:
                    normalized_tokens.append(token)
            if normalized_tokens:
                reactant_blocks.append(normalized_tokens)

        bucket["reactions"].append(
            {
                "id": legacy_reaction_id,
                "canonical_id": reaction.id,
                "name": reaction.name,
                "aliases": list(reaction.aliases),
                "description": reaction.description,
                "reactants": reactant_blocks,
                "required_roles": [
                    {
                        "role_id": role.role_id,
                        "required": bool(role.required),
                        "default_family_id": role.default_family_id,
                        "notes": role.notes,
                    }
                    for role in reaction.required_roles
                ],
                "metadata": copy.deepcopy(reaction.metadata),
                "source_ids": list(reaction.source_ids),
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
            canonical_id = reaction.get("canonical_id") or reaction_id

            metadata = copy.deepcopy(reaction)
            metadata["category_key"] = category_key
            metadata["category"] = category_label
            id_to_meta[reaction_id] = metadata

            register_alias(reaction_id, reaction_id)
            register_alias(canonical_id, reaction_id)
            register_alias(reaction.get("name", ""), reaction_id)
            for alias in reaction.get("aliases", []):
                register_alias(alias, reaction_id)

    for alias, canonical in CSV_REACTION_OVERRIDES.items():
        register_alias(alias, canonical)

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
    _, alias_to_id = _reaction_indices()
    return alias_to_id.get(label.strip().lower())


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


def required_reactant_categories(reaction_id: str) -> Optional[List[List[str]]]:
    """
    Return the reactant category requirements for a reaction id.

    The JSON schema stores reactant pairs as ``[[...electrophile...], [...nucleophile...]]``.
    """
    metadata = describe_reaction_type(reaction_id)
    if not metadata:
        return None
    reactants = metadata.get("reactants")
    if isinstance(reactants, list):
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
