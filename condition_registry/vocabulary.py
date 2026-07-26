"""Immutable public vocabulary for condition roles and families."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import FrozenSet

from .loader import load_taxonomy


@dataclass(frozen=True)
class ConditionDefinitionVocabulary:
    """Validated identifiers available to downstream declarative rules."""

    role_ids: FrozenSet[str]
    family_ids: FrozenSet[str]
    schema_version: str


@lru_cache(maxsize=1)
def load_condition_vocabulary() -> ConditionDefinitionVocabulary:
    """Return immutable role and family identifiers owned by the registry."""
    taxonomy = load_taxonomy()
    return ConditionDefinitionVocabulary(
        role_ids=frozenset(str(item["id"]) for item in taxonomy["roles"]),
        family_ids=frozenset(str(item["id"]) for item in taxonomy["families"]),
        schema_version=str(
            taxonomy.get("schema_version") or "roles_families.v1"
        ),
    )


__all__ = ["ConditionDefinitionVocabulary", "load_condition_vocabulary"]
