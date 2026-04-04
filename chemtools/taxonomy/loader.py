"""
Centralized taxonomy loading utilities.

Provides consistent, cached access to taxonomy data files across the codebase.
Single source of truth for all taxonomy loading operations.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, List

from .substituent_composer import load_organic_groups_with_compositions


# Taxonomy file paths (relative to this module)
_TAXONOMY_DIR = Path(__file__).resolve().parent / "data"
_CHEMTOOLS_DIR = Path(__file__).resolve().parents[1]
_RETRO_DATA_DIR = _CHEMTOOLS_DIR / "retro" / "data"
REACTION_TYPES_FILE = _TAXONOMY_DIR / "reaction_types.v4.0.json"
COMPOUND_LOGIC_FILE = _TAXONOMY_DIR / "compound_logic.json"
GROUP_LOGIC_FILE = _TAXONOMY_DIR / "group_logic.json"
ORGANIC_GROUPS_FILE = _TAXONOMY_DIR / "organic_groups.v1.3.json"
ORGANIC_COMPOUNDS_FILE = _TAXONOMY_DIR / "organic_compounds.v1.3.json"
COMPOUND_GENERATION_RULES_FILE = _TAXONOMY_DIR / "compound_generation_rules.v1.json"
COMPOUND_OVERRIDES_FILE = _TAXONOMY_DIR / "compound_overrides.v1.json"
SCAFFOLD_MOTIFS_FILE = _TAXONOMY_DIR / "scaffold_motifs.v1.3.json"
FEATURIZER_LOGIC_FILE = _TAXONOMY_DIR / "featurizer_logic.json"
SYNTHON_FILE = _TAXONOMY_DIR / "synthons.v1.json"
MOTIF_SCOPE_INDEX_FILE = _TAXONOMY_DIR / "motif_scope_index.v1.json"
RETRON_PATTERNS_FILE = _RETRO_DATA_DIR / "retron_patterns.json"
HTE_TEMPLATES_FILE = _RETRO_DATA_DIR / "hte_templates.json"
TRANSFORMATION_PATTERNS_FILE = _TAXONOMY_DIR / "transformation_patterns.json"


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
def load_organic_groups() -> Dict[str, Any]:
    """Load organic groups taxonomy JSON with composed substituent groups."""
    if not ORGANIC_GROUPS_FILE.exists():
        return {}

    try:
        return load_organic_groups_with_compositions(ORGANIC_GROUPS_FILE)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_organic_compounds() -> Dict[str, Any]:
    """Load organic compounds taxonomy JSON.
    
    Returns:
        Dictionary with compound type definitions.

    The canonical payload is generated from ``compound_generation_rules`` +
    ``compound_overrides``.
    """
    try:
        from .compound_catalog import build_documented_compound_catalog

        payload = build_documented_compound_catalog()
        return payload if isinstance(payload, dict) else {}
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_compound_generation_rules() -> Dict[str, Any]:
    """Load compound generation rules taxonomy JSON."""
    if not COMPOUND_GENERATION_RULES_FILE.exists():
        return {}
    try:
        with COMPOUND_GENERATION_RULES_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_compound_overrides() -> Dict[str, Any]:
    """Load explicit compound motif overrides taxonomy JSON."""
    if not COMPOUND_OVERRIDES_FILE.exists():
        return {}
    try:
        with COMPOUND_OVERRIDES_FILE.open("r", encoding="utf-8") as f:
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


@lru_cache(maxsize=1)
def load_featurizer_logic() -> Dict[str, Any]:
    """Load featurizer-specific chemistry logic JSON.

    Returns:
        Dictionary containing featurizer logic configuration.
    """
    if not FEATURIZER_LOGIC_FILE.exists():
        return {}

    try:
        with FEATURIZER_LOGIC_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_synthons() -> Dict[str, Any]:
    """Load synthon taxonomy JSON.

    Returns:
        Dictionary containing synthon class definitions.
    """
    if not SYNTHON_FILE.exists():
        return {}

    try:
        with SYNTHON_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_motif_scope_index() -> Dict[str, Any]:
    """Load prebuilt motif scope index JSON.

    Returns:
        Dictionary containing scope_map and scaffold parent metadata.
    """
    if not MOTIF_SCOPE_INDEX_FILE.exists():
        return {}

    try:
        with MOTIF_SCOPE_INDEX_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_retron_patterns() -> Dict[str, Any]:
    """Load retrosynthesis retron pattern library JSON."""
    if not RETRON_PATTERNS_FILE.exists():
        return {}

    try:
        with RETRON_PATTERNS_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_hte_templates() -> Dict[str, Any]:
    """Load HTE-backed retrosynthesis template library JSON."""
    if not HTE_TEMPLATES_FILE.exists():
        return {}

    try:
        with HTE_TEMPLATES_FILE.open("r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


@lru_cache(maxsize=1)
def load_transformation_patterns() -> Dict[str, Any]:
    """Load observable transformation/leaving-group pattern definitions JSON."""
    if not TRANSFORMATION_PATTERNS_FILE.exists():
        return {}

    try:
        with TRANSFORMATION_PATTERNS_FILE.open("r", encoding="utf-8") as f:
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
    load_organic_groups.cache_clear()
    load_organic_compounds.cache_clear()
    load_compound_generation_rules.cache_clear()
    load_compound_overrides.cache_clear()
    load_scaffold_motifs.cache_clear()
    load_featurizer_logic.cache_clear()
    load_synthons.cache_clear()
    load_motif_scope_index.cache_clear()
    load_retron_patterns.cache_clear()
    load_hte_templates.cache_clear()
    load_transformation_patterns.cache_clear()


__all__ = [
    "REACTION_TYPES_FILE",
    "COMPOUND_LOGIC_FILE",
    "GROUP_LOGIC_FILE",
    "ORGANIC_GROUPS_FILE",
    "ORGANIC_COMPOUNDS_FILE",
    "SCAFFOLD_MOTIFS_FILE",
    "FEATURIZER_LOGIC_FILE",
    "SYNTHON_FILE",
    "MOTIF_SCOPE_INDEX_FILE",
    "RETRON_PATTERNS_FILE",
    "HTE_TEMPLATES_FILE",
    "TRANSFORMATION_PATTERNS_FILE",
    "load_reaction_types_raw",
    "load_reaction_types_list",
    "load_reaction_types_dict",
    "load_compound_logic",
    "load_group_logic",
    "load_organic_groups",
    "load_organic_compounds",
    "load_compound_generation_rules",
    "load_compound_overrides",
    "load_scaffold_motifs",
    "load_featurizer_logic",
    "load_synthons",
    "load_motif_scope_index",
    "load_retron_patterns",
    "load_hte_templates",
    "load_transformation_patterns",
    "get_reaction_by_id",
    "get_reaction_metadata",
    "clear_taxonomy_cache",
]
