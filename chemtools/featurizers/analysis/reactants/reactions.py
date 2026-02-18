"""Reaction taxonomy accessors for reactant analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Set, Tuple

from ....taxonomy import loader, reaction_catalog


def get_reaction_type_definitions() -> Dict[str, Any]:
    """Load reaction type definitions keyed by canonical reaction ID."""
    return loader.load_reaction_types_dict()


def get_reaction_types_file() -> Path:
    """Return the canonical path to the reaction types taxonomy file."""
    return loader.REACTION_TYPES_FILE


def clear_reaction_type_cache() -> None:
    """Clear taxonomy caches used by reaction type lookups."""
    reaction_catalog.load_reaction_catalog.cache_clear()
    loader.clear_taxonomy_cache()


def list_reaction_type_ids() -> List[str]:
    """Return all canonical reaction type IDs."""
    return sorted(get_reaction_type_definitions().keys())


def normalize_reaction_type(reaction_id: str) -> Optional[str]:
    """Normalize an aliased reaction ID to its canonical taxonomy ID."""
    if not reaction_id:
        return None
    return reaction_catalog.resolve_reaction_type(reaction_id)


def describe_reaction_type(reaction_id: str) -> Optional[Dict[str, Any]]:
    """Return full metadata for the given reaction type ID."""
    rxn_id = normalize_reaction_type(reaction_id)
    if not rxn_id:
        return None
    return get_reaction_type_definitions().get(rxn_id)


def build_reaction_lookup() -> Tuple[Dict[str, str], Dict[str, Dict[str, Any]]]:
    """Build (alias -> canonical id, canonical id -> full data) lookups."""
    _, alias_map = reaction_catalog.load_reaction_catalog()
    return dict(alias_map), get_reaction_type_definitions()


def _iter_role_slots(rxn_data: Dict[str, Any]) -> Generator[Dict[str, Any], None, None]:
    reactants = rxn_data.get("reactants")
    if isinstance(reactants, dict):
        for slot in reactants.values():
            if isinstance(slot, dict):
                yield slot
            elif isinstance(slot, list):
                yield {"allowed": slot, "required": False}

    legacy = rxn_data.get("reactant_roles")
    if isinstance(legacy, dict):
        for slot in legacy.values():
            if isinstance(slot, dict):
                yield slot


def iter_reactions_for_category(
    category_id: str,
) -> Generator[Tuple[str, Dict[str, Any]], None, None]:
    """Iterate over all reactions that accept the given reactant category."""
    reactions = get_reaction_type_definitions()
    for rxn_id, rxn_data in reactions.items():
        if not isinstance(rxn_data, dict):
            continue
        for slot in _iter_role_slots(rxn_data):
            allowed = slot.get("allowed", [])
            if isinstance(allowed, list):
                for token in allowed:
                    if reaction_catalog.motif_tokens_compatible(category_id, str(token)):
                        yield rxn_id, rxn_data
                        break
                else:
                    continue
                break


def required_reactant_categories(reaction_id: str) -> Set[str]:
    """Return required reactant categories for a reaction definition."""
    rxn_id = normalize_reaction_type(reaction_id)
    if not rxn_id:
        return set()
    rxn_data = get_reaction_type_definitions().get(rxn_id)
    if not isinstance(rxn_data, dict):
        return set()

    categories: Set[str] = set()
    for slot in _iter_role_slots(rxn_data):
        required = bool(slot.get("required", False))
        if not required:
            continue
        allowed = slot.get("allowed", [])
        if isinstance(allowed, list):
            categories.update(str(item) for item in allowed if item)
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
