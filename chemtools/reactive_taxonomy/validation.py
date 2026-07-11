"""Validation for the isolated reactive-taxonomy data bundle."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List


DATA_DIR = Path(__file__).with_name("data")


def load_taxonomy_data() -> Dict[str, Any]:
    """Load all v1 taxonomy documents."""
    payload: Dict[str, Any] = {}
    for path in sorted(DATA_DIR.glob("*.v1.json")):
        with path.open("r", encoding="utf-8") as handle:
            payload[path.stem] = json.load(handle)
    return payload


def validate_taxonomy() -> List[str]:
    """Return validation errors; an empty list means the bundle is valid."""
    errors: List[str] = []
    try:
        payload = load_taxonomy_data()
    except (OSError, json.JSONDecodeError) as exc:
        return [f"taxonomy_load_failed:{exc}"]
    expected = {"contexts.v1", "handles.v1", "rendering.v1"}
    missing = expected - set(payload)
    if missing:
        errors.append(f"missing_taxonomy_files:{','.join(sorted(missing))}")
        return errors
    contexts = payload["contexts.v1"]
    tokens = contexts.get("tokens") or []
    precedence = contexts.get("precedence") or []
    if len(tokens) != len(set(tokens)):
        errors.append("duplicate_context_tokens")
    if set(tokens) != set(precedence):
        errors.append("context_precedence_mismatch")
    families = payload["handles.v1"].get("site_families") or {}
    required = {"leaving_group", "pronucleophile_XH", "transfer_group", "electrophilic_center"}
    if set(families) != required:
        errors.append("site_family_mismatch")
    return errors


__all__ = ["load_taxonomy_data", "validate_taxonomy"]
