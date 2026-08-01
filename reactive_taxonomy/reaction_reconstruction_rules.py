"""Validated declarative rules for grammar-free structural reconstruction."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Tuple


REACTION_RECONSTRUCTION_RULE_SCHEMA_VERSION = "1.0"
REACTION_RECONSTRUCTION_RULE_DEFINITION_VERSION = (
    "reaction_reconstruction_rules.v1"
)
_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_reconstruction_rules.v1.json"
)


@lru_cache(maxsize=1)
def load_reaction_reconstruction_rules() -> Tuple[Dict[str, Any], ...]:
    """Load structural site constraints and operator bindings."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != REACTION_RECONSTRUCTION_RULE_SCHEMA_VERSION:
        raise ValueError("Unsupported reaction-reconstruction rule schema")
    rules = payload.get("rules")
    if not isinstance(rules, list) or not rules:
        raise ValueError("Reaction-reconstruction rules must be a nonempty list")
    identifiers = [str(rule.get("id") or "") for rule in rules]
    if any(not identifier for identifier in identifiers) or len(identifiers) != len(
        set(identifiers)
    ):
        raise ValueError("Invalid or duplicate reaction-reconstruction rule ID")
    semantic_keys = {
        "compatible_named_families",
        "display_name",
        "generic_label",
        "named_family",
        "transformation_class",
    }
    for rule in rules:
        rule_id = str(rule["id"])
        if semantic_keys.intersection(rule):
            raise ValueError(
                f"Structural rule contains interpretation metadata: {rule_id}"
            )
        slots = rule.get("slots")
        bindings = rule.get("operator_slot_bindings")
        if not isinstance(slots, dict) or not slots:
            raise ValueError(f"Structural rule has no slots: {rule_id}")
        if not isinstance(bindings, dict) or not bindings:
            raise ValueError(f"Structural rule has no operator bindings: {rule_id}")
        if any(str(slot) not in slots for slot in bindings.values()):
            raise ValueError(f"Structural rule binds an unknown slot: {rule_id}")
        if not str(rule.get("operator_id") or ""):
            raise ValueError(f"Structural rule has no operator ID: {rule_id}")
    return tuple(
        sorted(
            (dict(rule) for rule in rules),
            key=lambda rule: (-int(rule.get("priority", 0)), str(rule["id"])),
        )
    )


__all__ = [
    "REACTION_RECONSTRUCTION_RULE_DEFINITION_VERSION",
    "REACTION_RECONSTRUCTION_RULE_SCHEMA_VERSION",
    "load_reaction_reconstruction_rules",
]
