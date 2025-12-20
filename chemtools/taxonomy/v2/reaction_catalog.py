"""
Reaction type catalog utilities for taxonomy v2.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


REACTION_TYPES_FILE = Path(__file__).resolve().parent / "reaction_types.v3.3.json"
_DEFAULT_SLOTS = ("electrophiles", "nucleophiles", "acids", "activators")


@dataclass(frozen=True)
class ReactionTypeDefinition:
    id: str
    name: str
    category: str
    aliases: List[str]
    description: Optional[str]
    reactants: Dict[str, List[str]]
    catalysts: List[str]
    conditions: Optional[str]
    metadata: Dict[str, Any]
    reference_reactions: List[str]
    notes: Optional[str]


def _load_payload(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Reaction types file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _clean_list(values: Any) -> List[str]:
    if not values:
        return []
    if isinstance(values, list):
        return [str(v) for v in values if isinstance(v, str) and v.strip()]
    return []


def _normalize_reactants(raw: Any) -> Dict[str, List[str]]:
    if not isinstance(raw, dict):
        return {}

    reactants: Dict[str, List[str]] = {}
    for slot in _DEFAULT_SLOTS:
        if slot in raw:
            reactants[slot] = _clean_list(raw.get(slot))

    for slot, values in raw.items():
        if slot in reactants:
            continue
        cleaned = _clean_list(values)
        if cleaned or isinstance(values, list):
            reactants[slot] = cleaned
    return reactants


@lru_cache(maxsize=1)
def load_reaction_catalog(
    path: Optional[Path] = None,
) -> Tuple[Dict[str, ReactionTypeDefinition], Dict[str, str]]:
    payload = _load_payload(path or REACTION_TYPES_FILE)
    reactions = payload.get("reaction_types") or []

    definitions: Dict[str, ReactionTypeDefinition] = {}
    alias_map: Dict[str, str] = {}

    def register_alias(alias: str, reaction_id: str) -> None:
        if not alias:
            return
        alias_map.setdefault(alias.strip().lower(), reaction_id)

    for entry in reactions:
        if not isinstance(entry, dict):
            continue
        rxn_id = str(entry.get("id") or "").strip()
        if not rxn_id:
            continue
        name = str(entry.get("name") or rxn_id)
        category = str(entry.get("category") or "")
        aliases = [str(a) for a in (entry.get("aliases") or []) if isinstance(a, str)]
        description = entry.get("description")
        reactants = _normalize_reactants(entry.get("reactants"))
        catalysts = [str(c) for c in (entry.get("catalysts") or []) if isinstance(c, str)]
        conditions = entry.get("conditions")
        metadata = dict(entry.get("metadata") or {})
        reference_reactions = [
            str(r) for r in (entry.get("reference_reactions") or []) if isinstance(r, str)
        ]
        notes = entry.get("notes")

        definitions[rxn_id] = ReactionTypeDefinition(
            id=rxn_id,
            name=name,
            category=category,
            aliases=aliases,
            description=description,
            reactants=reactants,
            catalysts=catalysts,
            conditions=conditions,
            metadata=metadata,
            reference_reactions=reference_reactions,
            notes=notes,
        )

        register_alias(rxn_id, rxn_id)
        register_alias(name, rxn_id)
        for alias in aliases:
            register_alias(alias, rxn_id)

    return definitions, alias_map


def list_reaction_type_ids() -> List[str]:
    definitions, _ = load_reaction_catalog()
    return sorted(definitions.keys())


def get_reaction_type(reaction_id: str) -> Optional[ReactionTypeDefinition]:
    if not reaction_id:
        return None
    definitions, _ = load_reaction_catalog()
    return definitions.get(reaction_id)


def resolve_reaction_type(label: Optional[str]) -> Optional[str]:
    if not label:
        return None
    _, alias_map = load_reaction_catalog()
    return alias_map.get(label.strip().lower())
