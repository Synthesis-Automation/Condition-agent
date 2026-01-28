"""Reaction type definitions and lookup utilities.

This module provides access to reaction type metadata from the taxonomy
(reaction_types.v4.0.json), including reaction type normalization, reaction
querying, and reactant requirement lookups.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Set, Tuple

from ....taxonomy import loader


def _load_reaction_types_raw() -> Dict[str, Any]:
    """Load raw reaction type definitions from reaction taxonomy.

    Returns a dictionary keyed by reaction type ID with full metadata.
    This is the single source of truth for reaction definitions.
    """
    return loader.load_reaction_types_dict()


def get_reaction_type_definitions() -> Dict[str, Any]:
    """Public API: Load reaction type definitions.

    Loads from reaction_types.v4.0.json in the taxonomy directory.
    """
    return _load_reaction_types_raw()


def get_reaction_types_file() -> Path:
    """Return the canonical path to the reaction types taxonomy file."""
    return loader.REACTION_TYPES_FILE


def clear_reaction_type_cache() -> None:
    """Clear all reaction type caches for testing or dynamic reloading."""
    _reaction_indices.cache_clear()


@lru_cache(maxsize=1)
def _reaction_indices() -> Dict[str, Any]:
    """Build lookup indices for reaction type queries (cached).

    Returns:
        Dictionary with keys:
        - 'by_id': reaction_id -> full reaction dict
        - 'by_alias': lowercase alias -> canonical reaction_id
        - 'by_family': family_id -> list of reaction_ids
        - 'families': family_id -> family metadata
    """
    reactions = _load_reaction_types_raw()

    by_id: Dict[str, dict] = {}
    by_alias: Dict[str, str] = {}
    by_family: Dict[str, List[str]] = {}
    families: Dict[str, dict] = {}

    for rxn_id, rxn_data in reactions.items():
        if not isinstance(rxn_data, dict):
            continue

        by_id[rxn_id] = rxn_data

        # Register canonical ID
        by_alias[rxn_id.lower()] = rxn_id

        # Register aliases
        for alias in rxn_data.get("aliases", []):
            if alias:
                by_alias[alias.lower()] = rxn_id

        # Track family membership
        family = rxn_data.get("family")
        if family:
            by_family.setdefault(family, []).append(rxn_id)
            if family not in families:
                # Extract family metadata from first member
                families[family] = {
                    "id": family,
                    "name": rxn_data.get("family_name", family),
                    "members": [],
                }
            families[family]["members"].append(rxn_id)

    return {
        "by_id": by_id,
        "by_alias": by_alias,
        "by_family": by_family,
        "families": families,
    }


def list_reaction_type_ids() -> List[str]:
    """Return all canonical reaction type IDs."""
    return list(_reaction_indices()["by_id"].keys())


def describe_reaction_type(reaction_id: str) -> Optional[Dict[str, Any]]:
    """Return full metadata for the given reaction type ID.

    Args:
        reaction_id: Canonical or aliased reaction type identifier

    Returns:
        Reaction metadata dict, or None if not found
    """
    reaction_id_normalized = normalize_reaction_type(reaction_id)
    if not reaction_id_normalized:
        return None
    return _reaction_indices()["by_id"].get(reaction_id_normalized)


def normalize_reaction_type(reaction_id: str) -> Optional[str]:
    """Normalize an aliased reaction ID to its canonical form.

    Args:
        reaction_id: Reaction identifier (canonical or alias)

    Returns:
        Canonical reaction type ID, or None if unrecognized
    """
    if not reaction_id:
        return None
    idx = _reaction_indices()
    # Try exact match first
    if reaction_id in idx["by_id"]:
        return reaction_id
    # Try case-insensitive alias lookup
    return idx["by_alias"].get(reaction_id.lower())


def build_reaction_lookup() -> Tuple[Dict[str, str], Dict[str, dict]]:
    """Build (alias -> canonical id, canonical id -> full data) lookups.

    Returns:
        Tuple of (alias_map, reaction_data) where:
        - alias_map: lowercase alias -> canonical reaction_id
        - reaction_data: canonical reaction_id -> full reaction dict
    """
    idx = _reaction_indices()
    return idx["by_alias"], idx["by_id"]


def iter_reactions_for_category(
    category_id: str,
) -> Generator[Tuple[str, Dict[str, Any]], None, None]:
    """Iterate over all reactions that accept the given reactant category.

    Args:
        category_id: Reactant category identifier (e.g., "organoboron")

    Yields:
        (reaction_id, reaction_metadata) tuples
    """
    reactions = _load_reaction_types_raw()
    for rxn_id, rxn_data in reactions.items():
        if not isinstance(rxn_data, dict):
            continue

        # Check if category appears in any role's allowed types
        for slot in rxn_data.get("reactant_roles", {}).values():
            if not isinstance(slot, dict):
                continue
            allowed = slot.get("allowed", [])
            if category_id in allowed:
                yield rxn_id, rxn_data
                break


def required_reactant_categories(reaction_id: str) -> Set[str]:
    """Return the set of required reactant categories for a reaction.

    Args:
        reaction_id: Canonical or aliased reaction type identifier

    Returns:
        Set of category IDs marked as required in the reaction definition
    """
    rxn_id = normalize_reaction_type(reaction_id)
    if not rxn_id:
        return set()

    rxn_data = _reaction_indices()["by_id"].get(rxn_id)
    if not rxn_data:
        return set()

    categories: Set[str] = set()
    for slot in rxn_data.get("reactant_roles", {}).values():
        if not isinstance(slot, dict):
            continue
        if slot.get("required", False):
            categories.update(slot.get("allowed", []))

    return categories


__all__ = [
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
