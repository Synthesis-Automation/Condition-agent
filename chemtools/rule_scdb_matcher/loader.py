from __future__ import annotations

"""Loader utilities for condition rule JSON payloads."""

import json
import re
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem

from .types import (
    CompiledSmarts,
    SchemeConditionDB,
    SchemeEntry,
    SelectorRule,
    SelectorRuleDB,
)


class RuleDatabaseError(RuntimeError):
    """Raised when a condition rule payload cannot be parsed."""


# Backwards compatibility with previous exception name
SchemeConditionDBError = RuleDatabaseError


def _compile_smarts(pattern: str) -> CompiledSmarts:
    mol = Chem.MolFromSmarts(pattern)
    if mol is None:
        raise RuleDatabaseError(f"Invalid SMARTS pattern: {pattern!r}")
    return CompiledSmarts(pattern=pattern, mol=mol, atom_count=mol.GetNumAtoms())


def _infer_entry_type(entry: dict[str, Any]) -> str:
    entry_type = entry.get("type")
    if entry_type:
        return str(entry_type)

    role = entry.get("role")
    if role == "default_condition":
        return "default_condition"

    return "scheme"


def _derive_reactant_smarts(entry: dict[str, Any]) -> tuple[str, ...]:
    smarts_list = entry.get("reactant_smarts") or []
    if smarts_list:
        return tuple(str(item) for item in smarts_list if item)

    rxn_smiles = entry.get("rxn_smiles_min") or entry.get("rxn_smiles")
    if isinstance(rxn_smiles, str) and ">>" in rxn_smiles:
        reactant_side = rxn_smiles.split(">>", 1)[0]
        candidates = [segment.strip() for segment in reactant_side.split(".") if segment.strip()]
        if candidates:
            return tuple(candidates)

    return tuple()


def _build_entry(entry: dict[str, Any]) -> SchemeEntry:
    entry_type = _infer_entry_type(entry)

    entry_id = entry.get("id")
    if not entry_id:
        raise RuleDatabaseError(f"Entry missing 'id': {entry}")

    reactant_smarts = _derive_reactant_smarts(entry)
    compiled = tuple(_compile_smarts(smarts) for smarts in reactant_smarts)

    applies_to = entry.get("applies_to")
    if applies_to is not None and not isinstance(applies_to, dict):
        raise RuleDatabaseError(
            f"Entry {entry_id!r} has unsupported 'applies_to' payload: {applies_to!r}"
        )

    env = entry.get("env") or {}
    feature_requirements = env.get("feature_requirements") if isinstance(env, dict) else None
    if feature_requirements is not None and not isinstance(feature_requirements, dict):
        raise RuleDatabaseError(
            f"Entry {entry_id!r} has unsupported 'feature_requirements': {feature_requirements!r}"
        )

    raw_copy = json.loads(json.dumps(entry))

    return SchemeEntry(
        id=str(entry_id),
        type=str(entry_type),
        name=entry.get("name"),
        reactant_smarts=reactant_smarts,
        applies_to=applies_to,
        feature_requirements=feature_requirements,
        priority=int(entry.get("priority", 0)),
        conditions=entry.get("conditions"),
        default_condition=entry.get("default_condition"),
        notes=tuple(entry.get("notes", []) or []),
        compiled_smarts=compiled,
        raw=raw_copy,
    )


def _normalize_reaction_type(payload: dict[str, Any]) -> str:
    reaction_type = payload.get("reaction_type")
    if isinstance(reaction_type, str) and reaction_type.strip():
        return reaction_type.strip()

    reaction = payload.get("reaction")
    if isinstance(reaction, str) and reaction.strip():
        slug = re.sub(r"[^A-Za-z0-9]+", "_", reaction).strip("_")
        return slug or reaction.strip()

    raise RuleDatabaseError("Database missing 'reaction_type'")


def _load_scheme_payload(payload: dict[str, Any], resolved: Path) -> SchemeConditionDB:
    schema_version = payload.get("schema_version")
    if not schema_version:
        raise RuleDatabaseError("Database missing 'schema_version'")

    reaction_type = _normalize_reaction_type(payload)
    updated_at = payload.get("updated_at")

    entries_data = payload.get("entries")
    if not isinstance(entries_data, list):
        raise RuleDatabaseError("Database 'entries' must be a list")

    entries = tuple(_build_entry(entry) for entry in entries_data)

    metadata = {
        key: value
        for key, value in payload.items()
        if key not in {"entries", "schema_version", "reaction_type", "reaction", "updated_at"}
    }

    return SchemeConditionDB(
        path=resolved,
        schema_version=str(schema_version),
        reaction_type=str(reaction_type),
        metadata=metadata,
        updated_at=str(updated_at) if updated_at is not None else None,
        entries=entries,
    )


def _ensure_dict(value: Any, context: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise RuleDatabaseError(f"Expected mapping for {context}, got {type(value).__name__}")
    return value


def _coerce_dict_sequence(items: Any, context: str) -> tuple[dict[str, Any], ...]:
    if not items:
        return tuple()
    if not isinstance(items, Iterable):
        raise RuleDatabaseError(f"Expected iterable for {context}, got {type(items).__name__}")
    result: list[dict[str, Any]] = []
    for index, item in enumerate(items):
        if not isinstance(item, dict):
            raise RuleDatabaseError(f"{context} entry {index} must be a mapping")
        result.append(dict(item))
    return tuple(result)


def _parse_rank_hint(selector: dict[str, Any], payload: dict[str, Any]) -> float | None:
    raw = payload.get("rank_hint", selector.get("rank_hint"))
    if raw is None:
        return None
    try:
        return float(raw)
    except (TypeError, ValueError) as exc:  # pragma: no cover - defensive branch
        raise RuleDatabaseError(f"Invalid rank_hint value for selector {selector.get('id')!r}: {raw!r}") from exc


def _build_selector(selector: dict[str, Any]) -> SelectorRule:
    selector_id = selector.get("id")
    if not selector_id:
        raise RuleDatabaseError(f"Selector missing 'id': {selector}")

    conditions = _ensure_dict(selector.get("if"), f"selector {selector_id!r} 'if'")
    payload = _ensure_dict(selector.get("then"), f"selector {selector_id!r} 'then'")
    rank_hint = _parse_rank_hint(selector, payload)

    raw_copy = json.loads(json.dumps(selector))

    return SelectorRule(
        id=str(selector_id),
        conditions=conditions,
        payload=payload,
        rank_hint=rank_hint,
        raw=raw_copy,
    )


def _load_selector_payload(payload: dict[str, Any], resolved: Path) -> SelectorRuleDB:
    schema_version = payload.get("schema_version")
    if not schema_version:
        raise RuleDatabaseError("Database missing 'schema_version'")

    reaction_type = _normalize_reaction_type(payload)

    selectors_data = payload.get("selectors")
    if not isinstance(selectors_data, list):
        raise RuleDatabaseError("Database 'selectors' must be a list")

    selectors = tuple(_build_selector(item) for item in selectors_data)

    defaults = _ensure_dict(payload.get("defaults"), "defaults")
    feature_schema = payload.get("feature_schema")
    if feature_schema is not None and not isinstance(feature_schema, dict):
        raise RuleDatabaseError("'feature_schema' must be a mapping when provided")

    guardrails = _coerce_dict_sequence(payload.get("guardrails"), "guardrails")
    priors = _coerce_dict_sequence(payload.get("priors"), "priors")
    repairs = _coerce_dict_sequence(payload.get("repairs"), "repairs")

    metadata = {
        key: value
        for key, value in payload.items()
        if key
        not in {
            "schema_version",
            "reaction_type",
            "reaction",
            "selectors",
            "defaults",
            "feature_schema",
            "guardrails",
            "priors",
            "repairs",
        }
    }

    return SelectorRuleDB(
        path=resolved,
        schema_version=str(schema_version),
        reaction_type=str(reaction_type),
        metadata=metadata,
        selectors=selectors,
        defaults=defaults,
        feature_schema=feature_schema,
        guardrails=guardrails,
        priors=priors,
        repairs=repairs,
    )


def load_db(path: str | Path) -> SchemeConditionDB | SelectorRuleDB:
    """Load and validate a condition rule JSON file."""

    resolved = Path(path).expanduser().resolve()
    with resolved.open("r", encoding="utf-8") as handle:
        payload: dict[str, Any] = json.load(handle)

    if "entries" in payload:
        return _load_scheme_payload(payload, resolved)
    if "selectors" in payload:
        return _load_selector_payload(payload, resolved)

    raise RuleDatabaseError("Unsupported database payload: expected 'entries' or 'selectors'")
