"""
Reaction type catalog utilities for taxonomy v2.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple


REACTION_TYPES_FILE = Path(__file__).resolve().parent / "data" / "reaction_types.v4.0.json"
COMPOUND_LOGIC_FILE = Path(__file__).resolve().parent / "data" / "compound_logic.json"
_DEFAULT_SLOTS = ("electrophiles", "nucleophiles", "acids", "activators", "substrate", "reagent")


@dataclass(frozen=True)
class ReactionTypeDefinition:
    id: str
    name: str
    category: str
    aliases: List[str]
    description: Optional[str]
    reactants: Dict[str, "SlotRequirement"]
    products: Dict[str, "SlotRequirement"]
    catalysts: List[str]
    conditions: Optional[str]
    metadata: Dict[str, Any]
    reference_reactions: List[str]
    notes: Optional[str]
    constraints: Dict[str, Any]


@dataclass(frozen=True)
class SlotRequirement:
    allowed: List[str]
    min_hits: int = 1
    min_reactants: int = 1

    def __contains__(self, item: object) -> bool:
        if not isinstance(item, str):
            return False
        return item in self.allowed


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


def _coerce_list(values: Any) -> List[str]:
    if isinstance(values, list):
        return [str(v) for v in values if isinstance(v, str) and v.strip()]
    if isinstance(values, str):
        value = values.strip()
        return [value] if value else []
    return []


def _dedupe(values: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    result: List[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def _expand_reactant_slot(
    values: Any,
    motif_sets: Mapping[str, Iterable[str]],
) -> SlotRequirement:
    min_hits = 1
    min_reactants = 1

    def resolve_value(v: str) -> List[str]:
        if v.startswith("@"):
            set_name = v[1:]
            entry = motif_sets.get(set_name)
            if entry:
                if isinstance(entry, dict):
                    return _coerce_list(entry.get("members"))
                return _coerce_list(entry)
        return [v]

    if isinstance(values, str):
        expanded = resolve_value(values)
        return SlotRequirement(allowed=_dedupe(expanded), min_hits=1, min_reactants=1)

    if isinstance(values, list):
        expanded = []
        for v in values:
            if isinstance(v, str):
                expanded.extend(resolve_value(v))
            else:
                expanded.append(str(v))
        return SlotRequirement(allowed=_dedupe(expanded), min_hits=1, min_reactants=1)

    if isinstance(values, dict):
        expanded: List[str] = []
        for set_name in _coerce_list(values.get("motif_sets") or values.get("motif_set")):
            entry = motif_sets.get(set_name)
            if isinstance(entry, dict):
                expanded.extend(_coerce_list(entry.get("members")))
            else:
                expanded.extend(_coerce_list(entry))

        for v in _coerce_list(values.get("include")):
            expanded.extend(resolve_value(v))

        exclude = set(_coerce_list(values.get("exclude")))
        if exclude:
            expanded = [value for value in expanded if value not in exclude]

        min_hits = int(values.get("min_hits") or 1)
        min_reactants = int(values.get("min_reactants") or 1)
        return SlotRequirement(
            allowed=_dedupe(expanded),
            min_hits=max(min_hits, 1),
            min_reactants=max(min_reactants, 1),
        )
    return SlotRequirement(allowed=[], min_hits=1, min_reactants=1)


def _normalize_reactants(
    raw: Any,
    motif_sets: Optional[Mapping[str, Iterable[str]]] = None,
) -> Dict[str, SlotRequirement]:
    if not isinstance(raw, dict):
        return {}

    motif_sets = motif_sets or {}
    reactants: Dict[str, SlotRequirement] = {}
    for slot in _DEFAULT_SLOTS:
        if slot in raw:
            reactants[slot] = _expand_reactant_slot(raw.get(slot), motif_sets)

    for slot, values in raw.items():
        if slot in reactants:
            continue
        cleaned = _expand_reactant_slot(values, motif_sets)
        if cleaned.allowed or isinstance(values, (list, dict)):
            reactants[slot] = cleaned
    return reactants


@lru_cache(maxsize=1)
def load_reaction_catalog(
    path: Optional[Path] = None,
) -> Tuple[Dict[str, ReactionTypeDefinition], Dict[str, str]]:
    payload = _load_payload(path or REACTION_TYPES_FILE)
    reactions = payload.get("reaction_types") or []
    motif_sets = payload.get("motif_sets") or {}

    # Load centralized compound logic if available
    logic_path = COMPOUND_LOGIC_FILE
    if logic_path.exists():
        try:
            with logic_path.open("r", encoding="utf-8") as h:
                logic_payload = json.load(h)
                logic_motif_sets = logic_payload.get("motif_sets") or {}
                # Merge logic_motif_sets into motif_sets (logic supplements)
                for k, v in logic_motif_sets.items():
                    if k not in motif_sets:
                        motif_sets[k] = v
        except Exception:
            pass

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
        reactants = _normalize_reactants(entry.get("reactants"), motif_sets)
        products = _normalize_reactants(entry.get("products"), motif_sets)
        catalysts = [str(c) for c in (entry.get("catalysts") or []) if isinstance(c, str)]
        conditions = entry.get("conditions")
        metadata = dict(entry.get("metadata") or {})
        constraints = dict(entry.get("constraints") or {})
        raw_references = entry.get("reference_reactions")
        if not raw_references:
            raw_references = entry.get("examples") or []
        reference_reactions = [str(r) for r in raw_references if isinstance(r, str)]
        notes = entry.get("notes")

        definitions[rxn_id] = ReactionTypeDefinition(
            id=rxn_id,
            name=name,
            category=category,
            aliases=aliases,
            description=description,
            reactants=reactants,
            products=products,
            catalysts=catalysts,
            conditions=conditions,
            metadata=metadata,
            reference_reactions=reference_reactions,
            notes=notes,
            constraints=constraints,
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
