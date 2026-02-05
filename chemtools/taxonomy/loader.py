"""
Centralized taxonomy loading utilities.

Provides consistent, cached access to taxonomy data files across the codebase.
Single source of truth for all taxonomy loading operations.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple


# Taxonomy file paths (relative to this module)
_TAXONOMY_DIR = Path(__file__).resolve().parent / "data"
REACTION_TYPES_FILE = _TAXONOMY_DIR / "reaction_types.v4.0.json"
COMPOUND_LOGIC_FILE = _TAXONOMY_DIR / "compound_logic.json"
GROUP_LOGIC_FILE = _TAXONOMY_DIR / "group_logic.json"
ORGANIC_COMPOUNDS_FILE = _TAXONOMY_DIR / "organic_compounds.v1.3.json"
SCAFFOLD_MOTIFS_FILE = _TAXONOMY_DIR / "scaffold_motifs.v1.3.json"


@lru_cache(maxsize=1)
def load_reaction_types_raw() -> Dict[str, Any]:
    """Load raw reaction types taxonomy JSON.
    
    Returns the complete JSON payload from reaction_types.v4.0.json
    including version, notes, and reaction_types list.
    
    Returns:
        Dictionary with keys: version, notes, reaction_types (list)
    """
    if not REACTION_TYPES_FILE.exists():
        raise FileNotFoundError(f"Reaction types file not found: {REACTION_TYPES_FILE}")
    
    with REACTION_TYPES_FILE.open("r", encoding="utf-8") as f:
        return json.load(f)


@lru_cache(maxsize=1)
def load_reaction_types_list() -> List[Dict[str, Any]]:
    """Load reaction types as a list.
    
    Returns:
        List of reaction type dictionaries from reaction_types.v4.0.json
    """
    payload = load_reaction_types_raw()
    return payload.get("reaction_types", [])


@lru_cache(maxsize=1)
def load_reaction_types_dict() -> Dict[str, Dict[str, Any]]:
    """Load reaction types as a dictionary keyed by reaction ID.
    
    Returns:
        Dictionary mapping reaction_id -> reaction definition
    """
    reaction_list = load_reaction_types_list()
    return {rxn["id"]: rxn for rxn in reaction_list if isinstance(rxn, dict) and "id" in rxn}


@lru_cache(maxsize=1)
def load_compound_logic() -> Dict[str, Any]:
    """Load compound logic taxonomy JSON.
    
    Returns:
        Dictionary with motif_sets and other compound logic configuration
    """
    if not COMPOUND_LOGIC_FILE.exists():
        return {}
    
    try:
        with COMPOUND_LOGIC_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_group_logic() -> Dict[str, Any]:
    """Load group logic taxonomy JSON.
    
    Returns:
        Dictionary with group_sets and optional group_elements configuration
    """
    if not GROUP_LOGIC_FILE.exists():
        return {}
    
    try:
        with GROUP_LOGIC_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_organic_compounds() -> Dict[str, Any]:
    """Load organic compounds taxonomy JSON.
    
    Returns:
        Dictionary with compound type definitions
    """
    if not ORGANIC_COMPOUNDS_FILE.exists():
        return {}
    
    try:
        with ORGANIC_COMPOUNDS_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_scaffold_motifs() -> Dict[str, Any]:
    """Load scaffold motifs taxonomy JSON.
    
    Returns:
        Dictionary with scaffold motif definitions
    """
    if not SCAFFOLD_MOTIFS_FILE.exists():
        return {}
    
    try:
        with SCAFFOLD_MOTIFS_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def get_reaction_by_id(reaction_id: str) -> Dict[str, Any] | None:
    """Get a single reaction definition by ID.
    
    Args:
        reaction_id: Reaction identifier (e.g., "Suzuki_miyaura")
        
    Returns:
        Reaction definition dict or None if not found
    """
    reactions = load_reaction_types_dict()
    return reactions.get(reaction_id)


def get_reaction_metadata(reaction_id: str, key: str | None = None) -> Any:
    """Get metadata for a reaction.
    
    Args:
        reaction_id: Reaction identifier
        key: Optional specific metadata key to retrieve
        
    Returns:
        Full metadata dict if key is None, otherwise specific value
    """
    reaction = get_reaction_by_id(reaction_id)
    if not reaction:
        return {} if key is None else None
    
    metadata = reaction.get("metadata", {})
    if key is None:
        return metadata
    return metadata.get(key)


def clear_taxonomy_cache() -> None:
    """Clear all cached taxonomy data.
    
    Use this when taxonomy files are modified and need to be reloaded.
    """
    load_reaction_types_raw.cache_clear()
    load_reaction_types_list.cache_clear()
    load_reaction_types_dict.cache_clear()
    load_compound_logic.cache_clear()
    load_group_logic.cache_clear()
    load_organic_compounds.cache_clear()
    load_scaffold_motifs.cache_clear()


__all__ = [
    "REACTION_TYPES_FILE",
    "COMPOUND_LOGIC_FILE",
    "GROUP_LOGIC_FILE",
    "ORGANIC_COMPOUNDS_FILE",
    "SCAFFOLD_MOTIFS_FILE",
    "load_reaction_types_raw",
    "load_reaction_types_list",
    "load_reaction_types_dict",
    "load_compound_logic",
    "load_group_logic",
    "load_organic_compounds",
    "load_scaffold_motifs",
    "get_reaction_by_id",
    "get_reaction_metadata",
    "clear_taxonomy_cache",
]
