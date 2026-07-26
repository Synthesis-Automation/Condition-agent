"""Immutable public vocabulary for condition roles and families."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, FrozenSet

from .loader import load_taxonomy

_DEFINITION_DIR = Path(__file__).with_name("definitions")


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


@lru_cache(maxsize=1)
def condition_registry_definition_versions() -> Dict[str, str]:
    """Return content-addressed versions for identity and role definitions."""
    versions = {}
    for path in sorted(_DEFINITION_DIR.glob("*"), key=lambda item: item.name):
        if not path.is_file() or path.name == "pending_substances.csv":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        versions[path.name] = f"sha256:{digest}"
    return versions


__all__ = [
    "ConditionDefinitionVocabulary",
    "condition_registry_definition_versions",
    "load_condition_vocabulary",
]
