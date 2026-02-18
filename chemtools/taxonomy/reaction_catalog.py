"""
Reaction type catalog utilities for taxonomy v2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple


REACTION_TYPES_FILE = Path(__file__).resolve().parent / "data" / "reaction_types.v4.0.json"
COMPOUND_LOGIC_FILE = Path(__file__).resolve().parent / "data" / "compound_logic.json"
SYNTHON_FILE = Path(__file__).resolve().parent / "data" / "synthons.v1.json"
_DEFAULT_SLOTS = ("electrophiles", "nucleophiles", "acids", "activators", "substrate", "reagent")
REACTION_CONSTRAINT_KEYS = (
    "include_reacted",
    "exclude_reacted",
    "include_formed",
    "exclude_formed",
    "include_bond_formed",
    "exclude_bond_formed",
    "include_bond_broken",
    "exclude_bond_broken",
    "min_reactant_slot_matches",
    "min_product_slot_matches",
)
REACTION_SYNTHON_SLOT_KEYS = (
    "synthon_sets",
    "synthon_set",
    "motif_sets",
    "motif_set",
    "include",
    "exclude",
    "min_hits",
    "min_reactants",
)


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
    redox_neutral: Optional[bool] = None
    synthons: Dict[str, "SlotRequirement"] = field(default_factory=dict)


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


def format_reaction_id_display(reaction_id: str) -> str:
    """Return a display-friendly string for a reaction id."""
    if not reaction_id:
        return ""
    tokens = reaction_id.split("_")
    if len(tokens) == 1:
        return reaction_id[:1].upper() + reaction_id[1:]
    formatted: List[str] = []
    for idx, token in enumerate(tokens):
        if not token:
            formatted.append(token)
            continue
        if len(token) == 1:
            formatted.append(token.upper())
            continue
        if token.isupper():
            formatted.append(token)
            continue
        if idx == 0:
            formatted.append(token[:1].upper() + token[1:])
        else:
            formatted.append(token)
    return "_".join(formatted)


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


@lru_cache(maxsize=1)
def _load_synthon_sets() -> Dict[str, List[str]]:
    if not SYNTHON_FILE.exists():
        return {}
    try:
        with SYNTHON_FILE.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}

    synthons = payload.get("synthons") or []
    by_role: Dict[str, List[str]] = {}
    all_ids: List[str] = []
    sets: Dict[str, List[str]] = {}

    for entry in synthons:
        if not isinstance(entry, dict):
            continue
        synthon_id = str(entry.get("id") or "").strip()
        role = str(entry.get("role") or "").strip().lower()
        if not synthon_id:
            continue
        sets[synthon_id] = [synthon_id]
        all_ids.append(synthon_id)
        if role:
            by_role.setdefault(role, []).append(synthon_id)

    for role, members in by_role.items():
        deduped = _dedupe(members)
        sets[role] = deduped
        sets[f"{role}_synthons"] = deduped
        if role == "electrophile":
            sets["electrophiles"] = deduped
        if role == "nucleophile":
            sets["nucleophiles"] = deduped

    sets["all"] = _dedupe(all_ids)

    custom_sets = payload.get("synthon_sets") or {}
    if isinstance(custom_sets, dict):
        for set_name, set_data in custom_sets.items():
            if not set_name:
                continue
            members = _coerce_list((set_data or {}).get("members")) if isinstance(set_data, dict) else _coerce_list(set_data)
            expanded: List[str] = []
            for member in members:
                if member.startswith("@"):
                    expanded.extend(sets.get(member[1:], []))
                else:
                    expanded.append(member)
            sets[str(set_name)] = _dedupe(expanded)

    return sets


def _expand_synthon_slot(
    values: Any,
    synthon_sets: Mapping[str, Iterable[str]],
) -> SlotRequirement:
    min_hits = 1
    min_reactants = 1

    def resolve_value(v: str) -> List[str]:
        if v.startswith("@"):
            return _coerce_list(synthon_sets.get(v[1:]))
        return [v]

    if isinstance(values, str):
        expanded = resolve_value(values)
        return SlotRequirement(allowed=_dedupe(expanded), min_hits=1, min_reactants=1)

    if isinstance(values, list):
        expanded = []
        for value in values:
            if isinstance(value, str):
                expanded.extend(resolve_value(value))
            else:
                expanded.append(str(value))
        return SlotRequirement(allowed=_dedupe(expanded), min_hits=1, min_reactants=1)

    if isinstance(values, dict):
        expanded: List[str] = []
        set_tokens = _coerce_list(values.get("synthon_sets") or values.get("synthon_set"))
        # Accept motif_set keys to keep schema evolution tolerant.
        if not set_tokens:
            set_tokens = _coerce_list(values.get("motif_sets") or values.get("motif_set"))
        for set_name in set_tokens:
            expanded.extend(_coerce_list(synthon_sets.get(set_name)))

        for value in _coerce_list(values.get("include")):
            expanded.extend(resolve_value(value))

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


def normalize_reaction_synthons(
    raw: Any,
    synthon_sets: Optional[Mapping[str, Iterable[str]]] = None,
) -> Dict[str, SlotRequirement]:
    if not isinstance(raw, dict):
        return {}
    synthon_sets = synthon_sets or {}
    normalized: Dict[str, SlotRequirement] = {}
    for slot, values in raw.items():
        cleaned = _expand_synthon_slot(values, synthon_sets)
        if cleaned.allowed or isinstance(values, (list, dict, str)):
            normalized[slot] = cleaned
    return normalized


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


def _to_non_negative_int(value: Any, default: int = 0) -> int:
    try:
        out = int(value)
    except Exception:
        out = int(default)
    return max(0, out)


def _coerce_optional_bool(value: Any) -> Optional[bool]:
    if isinstance(value, bool):
        return value
    return None


def normalize_reaction_constraints(raw: Any) -> Dict[str, Any]:
    """
    Normalize reaction constraints into a stable schema for all reaction families.

    Even if source JSON omits `constraints`, the loader provides this normalized
    payload with default values.
    """
    payload = dict(raw or {}) if isinstance(raw, dict) else {}
    return {
        "include_reacted": _dedupe(_coerce_list(payload.get("include_reacted"))),
        "exclude_reacted": _dedupe(_coerce_list(payload.get("exclude_reacted"))),
        "include_formed": _dedupe(_coerce_list(payload.get("include_formed"))),
        "exclude_formed": _dedupe(_coerce_list(payload.get("exclude_formed"))),
        "include_bond_formed": _dedupe(_coerce_list(payload.get("include_bond_formed"))),
        "exclude_bond_formed": _dedupe(_coerce_list(payload.get("exclude_bond_formed"))),
        "include_bond_broken": _dedupe(_coerce_list(payload.get("include_bond_broken"))),
        "exclude_bond_broken": _dedupe(_coerce_list(payload.get("exclude_bond_broken"))),
        "min_reactant_slot_matches": _to_non_negative_int(
            payload.get("min_reactant_slot_matches"),
            default=0,
        ),
        "min_product_slot_matches": _to_non_negative_int(
            payload.get("min_product_slot_matches"),
            default=0,
        ),
    }


@lru_cache(maxsize=1)
def _load_motif_scope_children() -> Dict[str, Set[str]]:
    """Load scope-map parent -> transitive children from taxonomy index."""
    try:
        from . import loader as taxonomy_loader

        payload = taxonomy_loader.load_motif_scope_index()
    except Exception:
        return {}

    raw_scope = payload.get("scope_map", {}) if isinstance(payload, dict) else {}
    if not isinstance(raw_scope, dict):
        return {}

    direct: Dict[str, Set[str]] = {}
    for key, value in raw_scope.items():
        parent = str(key).strip()
        if not parent:
            continue
        children = {
            str(item).strip()
            for item in (value or [])
            if str(item).strip() and str(item).strip() != parent
        }
        if children:
            direct[parent] = children

    if not direct:
        return {}

    transitive: Dict[str, Set[str]] = {}
    for parent in direct.keys():
        seen: Set[str] = set()
        stack = list(direct.get(parent, set()))
        while stack:
            child = stack.pop()
            if child in seen:
                continue
            seen.add(child)
            stack.extend(direct.get(child, set()))
        if seen:
            transitive[parent] = seen
    return transitive


@lru_cache(maxsize=1)
def _load_motif_scope_parents() -> Dict[str, Set[str]]:
    """Load inverse scope-map child -> transitive parents from taxonomy index."""
    children_map = _load_motif_scope_children()
    parents: Dict[str, Set[str]] = {}
    for parent, children in children_map.items():
        for child in children:
            parents.setdefault(child, set()).add(parent)
    return parents


def motif_tokens_compatible(left: str, right: str) -> bool:
    """
    Scope-aware token compatibility for taxonomy slot matching.

    A token is considered compatible if:
    - tokens are equal; or
    - either token is an in-scope parent/child of the other.
    """
    left_text = str(left or "").strip()
    right_text = str(right or "").strip()
    if not left_text or not right_text:
        return False
    if left_text == right_text:
        return True

    children_map = _load_motif_scope_children()
    if right_text in children_map.get(left_text, set()):
        return True
    if left_text in children_map.get(right_text, set()):
        return True

    parent_map = _load_motif_scope_parents()
    if right_text in parent_map.get(left_text, set()):
        return True
    if left_text in parent_map.get(right_text, set()):
        return True
    return False


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
    synthon_sets = _load_synthon_sets()

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
        raw_name = entry.get("name")
        name = str(raw_name) if raw_name else format_reaction_id_display(rxn_id)
        category = str(entry.get("category") or "")
        aliases = [str(a) for a in (entry.get("aliases") or []) if isinstance(a, str)]
        description = entry.get("description")
        reactants = _normalize_reactants(entry.get("reactants"), motif_sets)
        products = _normalize_reactants(entry.get("products"), motif_sets)
        catalysts = [str(c) for c in (entry.get("catalysts") or []) if isinstance(c, str)]
        conditions = entry.get("conditions")
        metadata = dict(entry.get("metadata") or {})
        constraints = normalize_reaction_constraints(entry.get("constraints"))
        redox_neutral = _coerce_optional_bool(entry.get("redox_neutral"))
        synthons = normalize_reaction_synthons(entry.get("synthons"), synthon_sets)
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
            redox_neutral=redox_neutral,
            synthons=synthons,
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
    definition = definitions.get(reaction_id)
    if definition:
        return definition
    resolved = resolve_reaction_type(reaction_id)
    if resolved:
        return definitions.get(resolved)
    return None


def resolve_reaction_type(label: Optional[str]) -> Optional[str]:
    if not label:
        return None
    _, alias_map = load_reaction_catalog()
    return alias_map.get(label.strip().lower())


__all__ = [
    "format_reaction_id_display",
    "load_reaction_catalog",
    "list_reaction_type_ids",
    "get_reaction_type",
    "resolve_reaction_type",
    "normalize_reaction_constraints",
    "normalize_reaction_synthons",
    "motif_tokens_compatible",
    "REACTION_CONSTRAINT_KEYS",
    "REACTION_SYNTHON_SLOT_KEYS",
]
