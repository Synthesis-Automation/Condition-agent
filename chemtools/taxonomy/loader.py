"""
Centralized taxonomy loading utilities.

Provides consistent, cached access to taxonomy data files across the codebase.
Single source of truth for all taxonomy loading operations.
"""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

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
            payload = json.load(f)
    except Exception:
        return {}
    return _expand_compound_logic_payload(payload)


def _dedupe_strs(values: Iterable[str]) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return out


def _compound_catalog_entries() -> List[Dict[str, str]]:
    payload = load_organic_compounds()
    compounds = payload.get("compounds", []) if isinstance(payload, dict) else []
    entries: List[Dict[str, str]] = []
    for entry in compounds:
        if not isinstance(entry, dict):
            continue
        compound_id = str(entry.get("id") or "").strip()
        group_a = str(entry.get("A") or "").strip()
        group_b = str(entry.get("B") or "").strip()
        if compound_id and group_a and group_b:
            entries.append({"id": compound_id, "A": group_a, "B": group_b})
    return entries


def _expand_group_set_members(
    set_id: str,
    group_sets: Mapping[str, Any],
    *,
    _seen: set[str] | None = None,
) -> List[str]:
    if not set_id or set_id not in group_sets:
        return []
    if _seen is None:
        _seen = set()
    if set_id in _seen:
        return []
    _seen.add(set_id)
    entry = group_sets.get(set_id) or {}
    members = entry.get("members", []) if isinstance(entry, dict) else []
    expanded: List[str] = []
    if not isinstance(members, list):
        return expanded
    for member in members:
        token = str(member or "").strip()
        if not token:
            continue
        if token in group_sets:
            expanded.extend(_expand_group_set_members(token, group_sets, _seen=_seen))
        else:
            expanded.append(token)
    return _dedupe_strs(expanded)


def _expand_compound_logic_axis(
    values: Any,
    *,
    group_sets: Mapping[str, Any],
) -> List[str]:
    if not values:
        return []
    items = values if isinstance(values, list) else [values]
    expanded: List[str] = []
    for item in items:
        if isinstance(item, str):
            token = item.strip()
            if token:
                expanded.append(token)
            continue
        if not isinstance(item, dict):
            continue
        set_id = str(item.get("set") or "").strip()
        if set_id:
            expanded.extend(_expand_group_set_members(set_id, group_sets))
            continue
        members = item.get("members")
        if isinstance(members, Sequence) and not isinstance(members, (str, bytes)):
            expanded.extend(str(member).strip() for member in members if str(member).strip())
    return _dedupe_strs(expanded)


def _expand_compound_logic_members(
    entry: Mapping[str, Any],
    *,
    group_sets: Mapping[str, Any],
    catalog_entries: Sequence[Mapping[str, str]],
) -> List[str]:
    members = _dedupe_strs(entry.get("members", []) or [])
    raw_rules = entry.get("member_rules", []) or []
    if not isinstance(raw_rules, list):
        raw_rules = []

    additions: List[str] = []
    removals: List[str] = []
    for rule in raw_rules:
        if not isinstance(rule, dict):
            continue
        additions.extend(_dedupe_strs(rule.get("members", []) or []))
        groups_a = _expand_compound_logic_axis(
            rule.get("A") or rule.get("groups_a"),
            group_sets=group_sets,
        )
        groups_b = _expand_compound_logic_axis(
            rule.get("B") or rule.get("groups_b"),
            group_sets=group_sets,
        )
        if groups_a and groups_b:
            allowed_a = set(groups_a)
            allowed_b = set(groups_b)
            for catalog_entry in catalog_entries:
                group_a = str(catalog_entry.get("A") or "").strip()
                group_b = str(catalog_entry.get("B") or "").strip()
                if group_a in allowed_a and group_b in allowed_b:
                    additions.append(str(catalog_entry.get("id") or "").strip())
        removals.extend(_dedupe_strs(rule.get("exclude", []) or []))

    expanded = _dedupe_strs([*members, *additions])
    if removals:
        removal_set = set(removals)
        expanded = [member for member in expanded if member not in removal_set]
    return expanded


def _expand_compound_logic_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(payload, dict):
        return {}
    motif_sets = payload.get("motif_sets") or {}
    if not isinstance(motif_sets, dict):
        return payload

    group_payload = load_group_logic()
    group_sets = group_payload.get("group_sets", {}) if isinstance(group_payload, dict) else {}
    catalog_entries = _compound_catalog_entries()

    expanded_payload = dict(payload)
    expanded_sets: Dict[str, Any] = {}
    for set_name, raw_entry in motif_sets.items():
        if not isinstance(raw_entry, dict):
            expanded_sets[str(set_name)] = raw_entry
            continue
        entry = dict(raw_entry)
        entry["members"] = _expand_compound_logic_members(
            entry,
            group_sets=group_sets if isinstance(group_sets, dict) else {},
            catalog_entries=catalog_entries,
        )
        expanded_sets[str(set_name)] = entry
    expanded_payload["motif_sets"] = expanded_sets
    return expanded_payload


@lru_cache(maxsize=1)
def load_compound_logic_sets() -> Dict[str, List[str]]:
    """Load expanded compound logic motif sets keyed by set id."""
    payload = load_compound_logic()
    motif_sets = payload.get("motif_sets", {}) if isinstance(payload, dict) else {}
    out: Dict[str, List[str]] = {}
    if not isinstance(motif_sets, dict):
        return out
    for set_name, entry in motif_sets.items():
        if isinstance(entry, dict):
            members = entry.get("members", []) or []
        elif isinstance(entry, list):
            members = entry
        else:
            members = []
        out[str(set_name)] = _dedupe_strs(members)
    return out


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

    The canonical payload is generated from ``compound_generation_rules``.
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
    load_compound_logic_sets.cache_clear()
    load_group_logic.cache_clear()
    load_organic_groups.cache_clear()
    load_organic_compounds.cache_clear()
    load_compound_generation_rules.cache_clear()
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
    "load_compound_logic_sets",
    "load_group_logic",
    "load_organic_groups",
    "load_organic_compounds",
    "load_compound_generation_rules",
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
