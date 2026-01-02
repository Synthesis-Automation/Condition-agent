"""
Loader for calculable feature specifications.

This version targets taxonomy v2 data. Reactant type features are generated
from ``organic_compounds.v1.3.json`` (with group templates) to avoid duplicating
SMARTS definitions. Optional property and derived overlays are loaded from
layered files when present.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

_BASE_FILE = "calculable_features.json"
_PROPERTIES_FILE = "calculable_features_properties.json"
_DERIVED_FILE = "calculable_features_derived.json"
_ORGANIC_COMPOUNDS_FILE = "organic_compounds.v1.3.json"
_ORGANIC_GROUPS_FILE = "organic_groups.v1.3.json"

_DEFAULT_TEMPLATES = {
    "single_bond": "{A}{B}",
    "via_oxygen": "{A}O{B}",
}


def _resolve_data_root(root: Optional[Path] = None) -> Path:
    if root is not None:
        return Path(root).resolve()
    return Path(__file__).resolve().parent / "data"


def _safe_load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object at {path}")
    return data


def _coerce_feature_list(value: Any, *, path: Path) -> list[dict[str, Any]]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"Expected 'features' to be a list at {path}")
    return [item for item in value if isinstance(item, dict)]


def _load_groups(path: Path) -> Dict[str, Dict[str, Any]]:
    payload = _safe_load_json(path)
    groups = payload.get("groups", [])
    return {
        group["id"]: group
        for group in groups
        if isinstance(group, dict) and group.get("id")
    }


def _load_compounds(path: Path) -> list[dict[str, Any]]:
    payload = _safe_load_json(path)
    compounds = payload.get("compounds") or []
    return [item for item in compounds if isinstance(item, dict)]


def _extract_compound_smarts(entry: dict[str, Any]) -> list[str]:
    smarts_any = entry.get("smarts_any") or entry.get("smarts")
    if isinstance(smarts_any, str):
        return [smarts_any]
    if isinstance(smarts_any, list):
        return [s for s in smarts_any if isinstance(s, str) and s.strip()]
    return []


def _build_reactant_spec(data_root: Path) -> Dict[str, Any]:
    compounds_path = data_root / _ORGANIC_COMPOUNDS_FILE
    groups_path = data_root / _ORGANIC_GROUPS_FILE
    if not compounds_path.exists() or not groups_path.exists():
        return {
            "version": "empty",
            "description": "Fallback empty calculable feature spec.",
            "features": [],
            "derived_shortcuts": [],
        }

    compounds = _load_compounds(compounds_path)
    groups = _load_groups(groups_path)
    templates = dict(_DEFAULT_TEMPLATES)
    features: list[dict[str, Any]] = []
    derived_shortcuts: list[dict[str, Any]] = []

    for entry in compounds:
        compound_id = entry.get("id") or ""
        if not compound_id:
            continue
        compound_name = entry.get("name") or compound_id
        compound_desc = entry.get("description") or ""
        group_b = entry.get("B") or ""

        smarts_list = _extract_compound_smarts(entry)
        if not smarts_list:
            template_id = entry.get("template") or ""
            template_format = templates.get(template_id)
            group_a = entry.get("A")
            group_b_id = entry.get("B")
            if template_format and group_a and group_b_id:
                group_a_record = groups.get(group_a, {})
                group_b_record = groups.get(group_b_id, {})
                a_smarts = group_a_record.get("smarts") or ""
                b_smarts = group_b_record.get("smarts") or ""
                if a_smarts and b_smarts:
                    smarts_list = [template_format.format(A=a_smarts, B=b_smarts)]

        if smarts_list:
            feature = {
                "token": f"{compound_id}_reactant",
                "type": "bool",
                "scope": "global",
                "category": "reactant_types",
                "detect": {"smarts_any": smarts_list},
                "why": compound_name,
                "metadata": {
                    "reactant_category": compound_id,
                    "reactant_member": compound_id,
                    "reactant_name": compound_name,
                    "category_name": compound_name,
                    "coupling_role": "other",
                    "legacy_taxonomy_id": compound_id,
                    "category_description": compound_desc,
                    "group": group_b,
                },
            }
            features.append(feature)

    return {
        "version": "organic_compounds_generated",
        "description": "Reactant feature spec generated from organic_compounds.v1.3.json.",
        "features": features,
        "derived_shortcuts": derived_shortcuts,
    }


@lru_cache(maxsize=4)
def load_calculable_feature_spec(root: Optional[Path] = None) -> Dict[str, Any]:
    data_root = _resolve_data_root(root)

    base_path = data_root / _BASE_FILE
    if base_path.exists():
        base = _safe_load_json(base_path)
        merged: Dict[str, Any] = dict(base)
        merged["features"] = _coerce_feature_list(base.get("features"), path=base_path)
        derived_shortcuts: list[dict[str, Any]] = []
        base_shortcuts = base.get("derived_shortcuts")
        if isinstance(base_shortcuts, list):
            derived_shortcuts = [item for item in base_shortcuts if isinstance(item, dict)]
    else:
        base = _build_reactant_spec(data_root)
        merged = dict(base)
        merged["features"] = _coerce_feature_list(base.get("features"), path=data_root / _ORGANIC_COMPOUNDS_FILE)
        derived_shortcuts = [
            item for item in (base.get("derived_shortcuts") or []) if isinstance(item, dict)
        ]

    for overlay_path in (data_root / _PROPERTIES_FILE, data_root / _DERIVED_FILE):
        if not overlay_path.exists():
            continue
        overlay = _safe_load_json(overlay_path)
        merged["features"].extend(_coerce_feature_list(overlay.get("features"), path=overlay_path))

        overlay_shortcuts = overlay.get("derived_shortcuts")
        if isinstance(overlay_shortcuts, list):
            derived_shortcuts.extend(item for item in overlay_shortcuts if isinstance(item, dict))

    merged["derived_shortcuts"] = derived_shortcuts

    if "changelog" in merged and merged["changelog"] is None:
        merged["changelog"] = []

    return merged


__all__ = ["load_calculable_feature_spec"]
