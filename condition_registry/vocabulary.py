"""Immutable public vocabulary for condition roles."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, FrozenSet

from .loader import load_role_definitions

_DEFINITION_DIR = Path(__file__).with_name("definitions")
_ACTIVE_DEFINITION_NAMES = {
    "recipe_templates.v1.json",
    "role_resolution.v2.json",
    "roles.v2.json",
    "substances.v2.jsonl",
    "synthesis_protocol.v1.schema.json",
}


@dataclass(frozen=True)
class ConditionDefinitionVocabulary:
    """Validated identifiers available to downstream declarative rules."""

    role_ids: FrozenSet[str]
    schema_version: str


@lru_cache(maxsize=1)
def load_condition_vocabulary() -> ConditionDefinitionVocabulary:
    """Return immutable role identifiers owned by the registry."""
    definitions = load_role_definitions()
    return ConditionDefinitionVocabulary(
        role_ids=frozenset(str(item["id"]) for item in definitions["roles"]),
        schema_version=str(definitions["schema_version"]),
    )


@lru_cache(maxsize=1)
def condition_registry_definition_versions() -> Dict[str, str]:
    """Return content-addressed versions for identity and role definitions."""
    versions = {}
    for path in sorted(_DEFINITION_DIR.glob("*"), key=lambda item: item.name):
        if not path.is_file() or path.name not in _ACTIVE_DEFINITION_NAMES:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        versions[path.name] = f"sha256:{digest}"
    return versions


__all__ = [
    "ConditionDefinitionVocabulary",
    "condition_registry_definition_versions",
    "load_condition_vocabulary",
]
