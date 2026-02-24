"""Reactant taxonomy loading and utilities.

This module handles loading and caching of reactant type definitions
from the organic_compounds taxonomy file.
"""

from __future__ import annotations

import copy
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from ....taxonomy import loader as taxonomy_loader
from .._registry import clear_registry_cache, get_registry

TAXONOMY_DATA_DIR = Path(__file__).resolve().parents[3] / "taxonomy" / "data"
REACTANT_TYPES_FILE = TAXONOMY_DATA_DIR / "organic_compounds.v1.3.json"


@lru_cache(maxsize=1)
def _load_reactant_types_raw() -> Dict[str, dict]:
    """Load reactant type definitions from registry or file."""
    registry = get_registry()
    if registry is None:
        try:
            payload = taxonomy_loader.load_organic_compounds()
        except Exception:
            payload = {}
        if isinstance(payload, dict) and payload.get("compounds"):
            return _load_reactant_types_from_compounds(payload.get("compounds") or [])
        return _load_reactant_types_from_file(REACTANT_TYPES_FILE)

    definitions: Dict[str, dict] = {}
    for reactant_id, reactant in registry.reactant_types.items():
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
        }
    return definitions


def _load_reactant_types_from_file(path: Path) -> Dict[str, dict]:
    """Load reactant types from JSON file."""
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
        }

    return definitions


def _load_reactant_types_from_compounds(compounds: Iterable[dict]) -> Dict[str, dict]:
    """Load reactant types from compounds-style JSON."""
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
                    "metadata": {"taxonomy_id": compound_id, "group": group},
                }
            ],
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
    """Build a lookup of aliases → canonical category IDs from taxonomy only."""
    definitions = _load_reactant_types_raw()
    alias_map: Dict[str, str] = {}
    for category, data in definitions.items():
        alias_map[category.lower()] = category
        for alias in data.get("aliases", []):
            alias_map[alias.lower()] = category
        for member in data.get("members", []):
            member_id = member.get("id", "")
            if member_id:
                alias_map[member_id.lower()] = category
            for alias in member.get("aliases", []):
                alias_map[alias.lower()] = category
    # Note: All aliases now sourced from taxonomy JSON files only
    return alias_map


def normalize_reactant_identifier(label: str) -> Optional[str]:
    """Return the canonical reactant category id for ``label``."""
    if not label:
        return None
    return _reactant_alias_index().get(label.strip().lower())


__all__ = [
    "REACTANT_TYPES_FILE",
    "get_reactant_type_definitions",
    "get_reactant_types_file",
    "clear_reactant_type_cache",
    "normalize_reactant_identifier",
    "_load_reactant_types_raw",
    "_reactant_alias_index",
]
