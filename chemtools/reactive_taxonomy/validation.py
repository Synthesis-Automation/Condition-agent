"""Validation for the isolated reactive-taxonomy data bundle."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from chemtools.core.smarts import compile_smarts


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
    patterns = payload["handles.v1"].get("patterns") or []
    pattern_ids = [str(pattern.get("id") or "") for pattern in patterns]
    if not patterns:
        errors.append("missing_handle_patterns")
    if len(pattern_ids) != len(set(pattern_ids)):
        errors.append("duplicate_handle_pattern_ids")
    for pattern in patterns:
        pattern_id = str(pattern.get("id") or "<missing>")
        if pattern.get("site_type") not in required:
            errors.append(f"invalid_pattern_site_type:{pattern_id}")
        smarts = str(pattern.get("smarts") or "")
        compiled = compile_smarts(smarts, validate=False)
        if compiled is None:
            errors.append(f"invalid_pattern_smarts:{pattern_id}")
            continue
        available_maps = {int(atom.GetAtomMapNum()) for atom in compiled.GetAtoms() if atom.GetAtomMapNum()}
        for role, raw_maps in (pattern.get("atom_roles") or {}).items():
            role_maps = raw_maps if isinstance(raw_maps, list) else [raw_maps]
            unknown_maps = {int(value) for value in role_maps} - available_maps
            if unknown_maps:
                errors.append(f"unknown_atom_map:{pattern_id}:{role}")
        for rule in pattern.get("suppresses") or []:
            if rule.get("site_type") not in required:
                errors.append(f"invalid_suppression_site_type:{pattern_id}")
            if rule.get("owned_role") not in (pattern.get("atom_roles") or {}):
                errors.append(f"invalid_suppression_role:{pattern_id}")
    return errors


__all__ = ["load_taxonomy_data", "validate_taxonomy"]
