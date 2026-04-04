"""Build canonical documented compound motifs from taxonomy generation rules.

This module is the migration path away from hand-maintained
``organic_compounds.v1.3.json``. It supports a rule-driven model:

- composable A/B motif entries from generation rules
- denylist filtering for unwanted combinations

The resulting payload intentionally matches the legacy compounds-style schema
(`{"compounds": [...]}`) so existing consumers can migrate incrementally.
"""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


_TAXONOMY_DIR = Path(__file__).resolve().parent / "data"
_GENERATION_RULES_FILE = _TAXONOMY_DIR / "compound_generation_rules.v1.json"
_GROUP_LOGIC_FILE = _TAXONOMY_DIR / "group_logic.json"


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _group_sets_index() -> Dict[str, Tuple[str, ...]]:
    payload = _load_json(_GROUP_LOGIC_FILE)
    out: Dict[str, Tuple[str, ...]] = {}
    raw_sets = payload.get("group_sets") or {}
    if not isinstance(raw_sets, dict):
        return out
    for set_id, entry in raw_sets.items():
        if not isinstance(entry, dict):
            continue
        members = entry.get("members") or []
        if not isinstance(members, list):
            continue
        cleaned = tuple(str(v).strip() for v in members if str(v).strip())
        if cleaned:
            out[str(set_id)] = cleaned
    return out


def _normalize_compound_id(entry: Mapping[str, Any]) -> str:
    explicit = str(entry.get("id") or "").strip()
    if explicit:
        return explicit
    group_a = str(entry.get("A") or "").strip()
    group_b = str(entry.get("B") or "").strip()
    if not group_a or not group_b:
        return ""
    return f"{group_a}{'' if group_b.startswith('-') else '-'}{group_b}"


def _normalize_template(entry: Mapping[str, Any], default_template: str) -> str:
    template = str(entry.get("template") or "").strip()
    return template or default_template


def _compound_identity_key(entry: Mapping[str, Any], *, default_template: str) -> Tuple[str, str, str]:
    group_a = str(entry.get("A") or "").strip()
    group_b = str(entry.get("B") or "").strip()
    template = _normalize_template(entry, default_template)
    return group_a, group_b, template


def _matches_deny_rule(
    entry: Mapping[str, Any],
    deny_rule: Mapping[str, Any],
    *,
    default_template: str,
) -> bool:
    deny_id = str(deny_rule.get("id") or "").strip()
    if deny_id and _normalize_compound_id(entry) == deny_id:
        return True

    deny_a = str(deny_rule.get("A") or "").strip()
    deny_b = str(deny_rule.get("B") or "").strip()
    deny_template = str(deny_rule.get("template") or "").strip()

    if deny_a and str(entry.get("A") or "").strip() != deny_a:
        return False
    if deny_b and str(entry.get("B") or "").strip() != deny_b:
        return False
    if deny_template and _normalize_template(entry, default_template) != deny_template:
        return False

    if deny_a or deny_b or deny_template:
        return True
    return False


def _clean_entry(entry: Mapping[str, Any], *, default_template: str, emit_ids: bool) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    compound_id = _normalize_compound_id(entry)
    if emit_ids and compound_id:
        out["id"] = compound_id

    # Preserve legacy-ish key order for readability.
    template = _normalize_template(entry, default_template)
    if "smarts" not in entry and "smarts_any" not in entry:
        if template:
            out["template"] = template
    else:
        if template and str(entry.get("template") or "").strip():
            out["template"] = template

    for key in ("A", "B", "smarts", "smarts_any", "description", "priority", "reactivity_weight", "aliases", "anchors"):
        if key not in entry:
            continue
        value = entry.get(key)
        if value is None:
            continue
        out[key] = deepcopy(value)

    # Ensure A/B survive when template omitted but ID is generated from them.
    if "A" not in out and entry.get("A") is not None:
        out["A"] = deepcopy(entry.get("A"))
    if "B" not in out and entry.get("B") is not None:
        out["B"] = deepcopy(entry.get("B"))

    return out


def _expand_member_list(values: Any, *, group_sets: Mapping[str, Tuple[str, ...]]) -> List[str]:
    if not values:
        return []
    if isinstance(values, str):
        token = values.strip()
        return [token] if token else []
    if not isinstance(values, list):
        return []

    out: List[str] = []
    for item in values:
        if isinstance(item, str):
            token = item.strip()
            if token:
                out.append(token)
            continue
        if not isinstance(item, dict):
            continue
        set_id = str(item.get("set") or "").strip()
        if set_id and set_id in group_sets:
            out.extend(group_sets[set_id])
            continue
        # Alternate explicit syntax: {"members": [...]}
        members = item.get("members")
        if isinstance(members, list):
            out.extend(str(v).strip() for v in members if str(v).strip())

    seen: Set[str] = set()
    deduped: List[str] = []
    for token in out:
        if token in seen:
            continue
        seen.add(token)
        deduped.append(token)
    return deduped


def _expand_rule_pairs(
    rule: Mapping[str, Any],
    *,
    default_template: str,
    group_sets: Mapping[str, Tuple[str, ...]],
) -> List[Dict[str, Any]]:
    template = str(rule.get("template") or "").strip() or default_template
    scaffolds = _expand_member_list(rule.get("scaffolds"), group_sets=group_sets)
    substituents = _expand_member_list(rule.get("substituents"), group_sets=group_sets)
    if not scaffolds or not substituents:
        return []

    deny_pairs = rule.get("deny_pairs") or []
    if not isinstance(deny_pairs, list):
        deny_pairs = []

    base_meta = rule.get("metadata") or {}
    if not isinstance(base_meta, dict):
        base_meta = {}

    entries: List[Dict[str, Any]] = []
    for group_a in scaffolds:
        for group_b in substituents:
            entry: Dict[str, Any] = {"A": group_a, "B": group_b, "template": template}
            entry.update(deepcopy(base_meta))
            if any(
                _matches_deny_rule(entry, deny_rule, default_template=default_template)
                for deny_rule in deny_pairs
                if isinstance(deny_rule, dict)
            ):
                continue
            entries.append(entry)
    return entries


def _build_from_rules(
    generation_rules: Mapping[str, Any],
    *,
    emit_ids: bool,
) -> Dict[str, Any]:
    defaults = generation_rules.get("defaults") or {}
    if not isinstance(defaults, dict):
        defaults = {}
    default_template = str(defaults.get("template") or "single_bond").strip() or "single_bond"
    group_sets = _group_sets_index()

    generated_entries: List[Dict[str, Any]] = []

    include_pairs = generation_rules.get("include_pairs") or []
    if isinstance(include_pairs, list):
        for entry in include_pairs:
            if isinstance(entry, dict):
                generated_entries.append(dict(entry))

    include_rules = generation_rules.get("include_rules") or []
    if isinstance(include_rules, list):
        for rule in include_rules:
            if not isinstance(rule, dict):
                continue
            generated_entries.extend(
                _expand_rule_pairs(rule, default_template=default_template, group_sets=group_sets)
            )

    deny_rules = generation_rules.get("deny_pairs") or []
    if not isinstance(deny_rules, list):
        deny_rules = []

    filtered_entries: List[Dict[str, Any]] = []
    seen_keys: Set[Tuple[str, str, str]] = set()
    seen_ids: Set[str] = set()
    for entry in generated_entries:
        if any(
            _matches_deny_rule(entry, deny_rule, default_template=default_template)
            for deny_rule in deny_rules
            if isinstance(deny_rule, dict)
        ):
            continue

        compound_id = _normalize_compound_id(entry)
        ident_key = _compound_identity_key(entry, default_template=default_template)
        if compound_id:
            if compound_id in seen_ids:
                continue
            seen_ids.add(compound_id)
        elif ident_key in seen_keys:
            continue
        seen_keys.add(ident_key)
        filtered_entries.append(_clean_entry(entry, default_template=default_template, emit_ids=emit_ids))

    return {
        "version": str(generation_rules.get("legacy_version") or "v1.3"),
        "naming_set": str(generation_rules.get("naming_set") or "chemist_style"),
        "notes": str(
            generation_rules.get("notes")
            or "Generated documented compound motifs from compound_generation_rules."
        ),
        "compounds": filtered_entries,
    }


def has_compound_catalog_rules() -> bool:
    """Return True when the split catalog files are present."""
    return _GENERATION_RULES_FILE.exists()


def build_documented_compound_catalog(
    *,
    emit_ids: bool = True,
) -> Dict[str, Any]:
    """Return the canonical documented compound motif catalog payload.

    Args:
        emit_ids: Materialize ``id`` for all A/B motifs to simplify downstream
            consumers that previously depended on explicit IDs.
    """
    generation_rules = _load_json(_GENERATION_RULES_FILE)
    if generation_rules:
        return _build_from_rules(generation_rules, emit_ids=emit_ids)
    return {}


__all__ = [
    "build_documented_compound_catalog",
    "has_compound_catalog_rules",
]
