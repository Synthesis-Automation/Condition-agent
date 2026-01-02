"""
Loader for calculable feature specifications.

This version targets taxonomy v2 data. Reactant type features are generated
from ``reactant_types.json`` to avoid duplicating SMARTS definitions. Optional
property and derived overlays are loaded from layered files when present.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

_BASE_FILE = "calculable_features.json"
_PROPERTIES_FILE = "calculable_features_properties.json"
_DERIVED_FILE = "calculable_features_derived.json"
_REACTANT_TYPES_FILE = "reactant_types.json"

_ROLE_MAPPING = {
    "ArX*": "electrophile",
    "HetAr-X": "electrophile",
    "VinylX*": "electrophile",
    "Alkyl-X": "electrophile",
    "Allylic-X": "electrophile",
    "Benzylic-X": "electrophile",
    "RSO2Cl": "electrophile",
    "Acyl-electrophile": "electrophile",
    "ArB*": "nucleophile",
    "RB*": "nucleophile",
    "R-M": "nucleophile",
    "RMgX": "nucleophile",
    "RZnX": "nucleophile",
    "RLi": "nucleophile",
    "RNH2/R2NH": "nucleophile",
    "ArNH2/Ar2NH": "nucleophile",
    "ROH": "nucleophile",
    "ArOH": "nucleophile",
    "RSH": "nucleophile",
    "Alkyne": "both",
    "Alkene": "both",
    "Aldehyde": "electrophile",
    "Ketone": "electrophile",
    "Ester": "electrophile",
    "RCO2H": "nucleophile",
    "Amide": "nucleophile",
    "Metal-H": "reductant",
    "H2": "reductant",
    "Oxidant": "oxidant",
    "Diene": "dienophile",
    "Dienophile": "dienophile",
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


def _load_reactant_types(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, list):
        raise ValueError(f"Expected a JSON list at {path}")
    return [item for item in payload if isinstance(item, dict)]


def _build_reactant_spec(data_root: Path) -> Dict[str, Any]:
    reactant_path = data_root / _REACTANT_TYPES_FILE
    if not reactant_path.exists():
        return {
            "version": "empty",
            "description": "Fallback empty calculable feature spec.",
            "features": [],
            "derived_shortcuts": [],
        }

    reactant_data = _load_reactant_types(reactant_path)
    features: list[dict[str, Any]] = []
    derived_shortcuts: list[dict[str, Any]] = []

    for category in reactant_data:
        cat_id = category.get("id") or ""
        cat_name = category.get("name") or cat_id
        cat_desc = category.get("description") or ""
        members = category.get("members") or []

        member_tokens: list[str] = []
        for member in members:
            member_id = member.get("id") or ""
            smarts = member.get("smarts") or ""
            if not member_id or not smarts:
                continue
            member_name = member.get("name") or member_id
            meta = member.get("metadata") or {}
            token = meta.get("feature_token") or f"{member_id}_reactant"
            legacy_id = meta.get("legacy_taxonomy_id") or member_id
            coupling_role = meta.get("coupling_role") or _ROLE_MAPPING.get(cat_id, "other")
            feature = {
                "token": token,
                "type": "bool",
                "scope": "global",
                "category": "reactant_types",
                "detect": {"smarts_any": [smarts]},
                "why": f"{member_name} - {cat_name}",
                "metadata": {
                    "reactant_category": cat_id,
                    "reactant_member": member_id,
                    "reactant_name": member_name,
                    "category_name": cat_name,
                    "coupling_role": coupling_role,
                    "legacy_taxonomy_id": legacy_id,
                    "category_description": cat_desc,
                },
            }
            features.append(feature)
            member_tokens.append(token)

        if member_tokens:
            safe_cat_id = (
                cat_id.replace("*", "")
                .replace("-", "_")
                .replace("/", "_")
            )
            derived_shortcuts.append(
                {
                    "token": f"{safe_cat_id}_reactant",
                    "derive": " OR ".join(member_tokens),
                    "why": f"Category-level: {cat_name}",
                    "metadata": {
                        "reactant_category": cat_id,
                        "category_name": cat_name,
                        "is_category_level": True,
                        "member_count": len(member_tokens),
                    },
                }
            )

    return {
        "version": "reactant_types_generated",
        "description": "Reactant feature spec generated from reactant_types.json.",
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
        merged["features"] = _coerce_feature_list(base.get("features"), path=data_root / _REACTANT_TYPES_FILE)
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
