"""
Reactant and reaction type utilities shared across ChemTools and HTE pipelines.

This module centralises the SMARTS-driven reactant classification helpers and the
reaction type taxonomies that previously lived under ``data-processor/HTE_data``.
It provides:

* Cached loaders for the JSON definitions.
* Reactant SMARTS classification that degrades gracefully when RDKit is absent.
* Lookup/normalisation helpers for datasets (e.g., z-score CSV processing).
* Reaction type alias normalisation and metadata introspection.
"""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from chemtools.util import rdkit_helpers

DATA_DIR = Path(__file__).resolve().parent / "data"
REACTANT_TYPES_FILE = DATA_DIR / "reactant_types.json"
REACTION_TYPES_FILE = DATA_DIR / "reaction_types.json"

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


@lru_cache(maxsize=1)
def _load_reactant_types_raw() -> Dict[str, dict]:
    with REACTANT_TYPES_FILE.open("r", encoding="utf-8") as fh:
        return json.load(fh)


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


@lru_cache(maxsize=1)
def _reactant_alias_index() -> Tuple[Dict[str, str], Dict[str, str]]:
    """Build (alias -> id, id -> category) mappings for reactant types."""
    definitions = _load_reactant_types_raw()
    alias_to_id: Dict[str, str] = {}
    id_to_category: Dict[str, str] = {}

    def register(alias: str, canonical: str) -> None:
        if not alias:
            return
        alias_to_id[alias.lower()] = canonical

    for category_key, category_data in definitions.items():
        register(category_key, category_key)
        id_to_category[category_key] = category_key

        # Category description can act as an alias when unique.
        description = category_data.get("description")
        if isinstance(description, str) and description:
            register(description, category_key)

        for member in category_data.get("members", []):
            member_id = member.get("id")
            if not member_id:
                continue
            register(member_id, member_id)
            id_to_category[member_id] = category_key

            # Names are useful aliases (case-insensitive).
            name = member.get("name")
            if isinstance(name, str) and name:
                register(name, member_id)

    for alias, canonical in CSV_REACTANT_OVERRIDES.items():
        register(alias, canonical)

    return alias_to_id, id_to_category


def build_reactant_lookup() -> Tuple[Dict[str, str], Dict[str, str]]:
    """Expose the cached lookup tables (alias -> id, id -> category)."""
    alias_to_id, id_to_category = _reactant_alias_index()
    return copy.deepcopy(alias_to_id), copy.deepcopy(id_to_category)


def normalize_reactant_identifier(value: str) -> Optional[str]:
    """Normalise a raw reactant label to a known member/category id."""
    if not value:
        return None
    alias_to_id, _ = _reactant_alias_index()
    return alias_to_id.get(value.strip().lower())


def reactant_category_for(identifier: str) -> Optional[str]:
    """Resolve the category for a given member/category identifier."""
    if not identifier:
        return None
    _, id_to_category = _reactant_alias_index()
    return id_to_category.get(identifier)


def iter_reactant_matches(
    smiles: str, reactant_types: Optional[Dict[str, dict]] = None
) -> List[ReactantMatch]:
    """
    Return all SMARTS matches for the provided SMILES string.

    Returns an empty list when RDKit is unavailable or parsing fails.
    """
    mol = rdkit_helpers.parse_smiles(smiles)
    if mol is None:
        return []

    definitions = reactant_types if reactant_types is not None else _load_reactant_types_raw()
    matches: List[ReactantMatch] = []

    for category, data in definitions.items():
        category_smarts = data.get("smarts")
        category_description = data.get("description", "")
        group = data.get("group", "")
        is_general = category in GENERAL_REACTANT_CATEGORIES

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
                        category=category,
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
    with REACTION_TYPES_FILE.open("r", encoding="utf-8") as fh:
        return json.load(fh)


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
