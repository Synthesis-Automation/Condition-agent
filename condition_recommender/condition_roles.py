"""Translate reaction interpretation into generic condition-role preferences."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Tuple

from condition_registry import load_condition_vocabulary

_RULES_PATH = Path(__file__).with_name("definitions") / "condition_role_context.v1.json"


@lru_cache(maxsize=1)
def load_condition_role_context() -> Dict[str, Any]:
    """Load and validate recommender-owned reaction-to-role preferences."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.0":
        raise ValueError("Unsupported condition role context schema")
    known_roles = set(load_condition_vocabulary().role_ids)
    referenced = {
        str(role)
        for roles in rules.get("transformation_role_preferences", {}).values()
        for role in roles
    }
    unknown = referenced - known_roles
    if unknown:
        raise ValueError(f"Condition role context has unknown roles: {sorted(unknown)}")
    return rules


def preferred_roles_for_reaction(transformation_class: str | None) -> Tuple[str, ...]:
    """Return generic role preferences for one interpreted transformation."""
    return tuple(
        str(value)
        for value in load_condition_role_context()
        .get("transformation_role_preferences", {})
        .get(transformation_class or "", ())
    )


__all__ = ["load_condition_role_context", "preferred_roles_for_reaction"]
