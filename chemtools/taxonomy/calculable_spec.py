"""
Loader for calculable feature specifications.

This version targets taxonomy v2 data. When the layered calculable feature
files are missing, it falls back to the generated reactant feature list
in scripts/reactant_features_generated.json so reactant classification
still works without archived taxonomy assets.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

_BASE_FILE = "calculable_features.json"
_PROPERTIES_FILE = "calculable_features_properties.json"
_DERIVED_FILE = "calculable_features_derived.json"


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


def _fallback_reactant_spec() -> Dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    fallback_path = root / "scripts" / "reactant_features_generated.json"
    if not fallback_path.exists():
        return {
            "version": "empty",
            "description": "Fallback empty calculable feature spec.",
            "features": [],
            "derived_shortcuts": [],
        }
    try:
        payload = _safe_load_json(fallback_path)
    except Exception:
        return {
            "version": "empty",
            "description": "Fallback empty calculable feature spec.",
            "features": [],
            "derived_shortcuts": [],
        }
    features = _coerce_feature_list(payload.get("features"), path=fallback_path)
    derived_shortcuts = payload.get("derived_shortcuts")
    if not isinstance(derived_shortcuts, list):
        derived_shortcuts = []
    derived_shortcuts = [item for item in derived_shortcuts if isinstance(item, dict)]
    return {
        "version": "reactant_features_generated",
        "description": "Reactant feature fallback spec.",
        "features": features,
        "derived_shortcuts": derived_shortcuts,
    }


@lru_cache(maxsize=4)
def load_calculable_feature_spec(root: Optional[Path] = None) -> Dict[str, Any]:
    data_root = _resolve_data_root(root)

    base_path = data_root / _BASE_FILE
    if not base_path.exists():
        return _fallback_reactant_spec()

    base = _safe_load_json(base_path)
    merged: Dict[str, Any] = dict(base)
    merged["features"] = _coerce_feature_list(base.get("features"), path=base_path)

    derived_shortcuts: list[dict[str, Any]] = []
    base_shortcuts = base.get("derived_shortcuts")
    if isinstance(base_shortcuts, list):
        derived_shortcuts = [item for item in base_shortcuts if isinstance(item, dict)]

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
