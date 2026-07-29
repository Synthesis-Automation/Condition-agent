"""Validated loading for versioned reactivity descriptor definitions."""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Tuple


_DEFINITIONS_DIR = Path(__file__).resolve().parent.parent / "definitions"
_DEFINITION_FILES = (
    "aromatic_systems.v1.json",
    "reactivity_descriptor_rules.v1.json",
    "reactivity_rendering.v1.json",
)


@lru_cache(maxsize=1)
def load_reactivity_descriptor_definitions() -> Dict[str, Dict[str, Any]]:
    """Load and minimally validate descriptor definitions."""
    payload: Dict[str, Dict[str, Any]] = {}
    for filename in _DEFINITION_FILES:
        path = _DEFINITIONS_DIR / filename
        with path.open("r", encoding="utf-8-sig") as handle:
            value = json.load(handle)
        if value.get("schema_version") != "1.0":
            raise ValueError(f"unsupported descriptor definition schema:{filename}")
        if not value.get("definition_id"):
            raise ValueError(f"missing descriptor definition id:{filename}")
        payload[filename] = dict(value)
    rules = payload["reactivity_descriptor_rules.v1.json"]
    steric = rules.get("steric") or {}
    electronic = rules.get("electronic") or {}
    if int(steric.get("radius", 0)) < 1:
        raise ValueError("invalid descriptor steric radius")
    if int(electronic.get("radius", 0)) < 1:
        raise ValueError("invalid descriptor electronic radius")
    for key in ("bins", "ortho_burden_bins"):
        values = steric.get(key) or {}
        if not values or any(not 0.0 <= float(value) <= 1.0 for value in values.values()):
            raise ValueError(f"invalid steric descriptor bins:{key}")
    if not electronic.get("activation_axes"):
        raise ValueError("missing electronic activation axes")
    return payload


@lru_cache(maxsize=1)
def descriptor_definition_versions() -> Tuple[Tuple[str, str], ...]:
    """Return deterministic schema and content hashes for profile provenance."""
    versions = []
    for filename in _DEFINITION_FILES:
        path = _DEFINITIONS_DIR / filename
        raw = path.read_bytes()
        with path.open("r", encoding="utf-8-sig") as handle:
            payload = json.load(handle)
        versions.append(
            (
                filename,
                f"{payload.get('schema_version', 'unknown')}@sha256:"
                f"{hashlib.sha256(raw).hexdigest()[:16]}",
            )
        )
    return tuple(sorted(versions))


def descriptor_rules() -> Dict[str, Any]:
    """Return the validated core descriptor-rule document."""
    return dict(
        load_reactivity_descriptor_definitions()[
            "reactivity_descriptor_rules.v1.json"
        ]
    )


def aromatic_system_rules() -> Dict[str, Any]:
    """Return validated aromatic-system classification rules."""
    return dict(
        load_reactivity_descriptor_definitions()["aromatic_systems.v1.json"]
    )


def rendering_rules() -> Dict[str, Any]:
    """Return validated chemist-facing descriptor rendering rules."""
    return dict(
        load_reactivity_descriptor_definitions()["reactivity_rendering.v1.json"]
    )


__all__ = [
    "aromatic_system_rules",
    "descriptor_definition_versions",
    "descriptor_rules",
    "load_reactivity_descriptor_definitions",
    "rendering_rules",
]
