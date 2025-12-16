"""
Layered loader for the calculable feature specification.

The feature system is intentionally split into layers (see
`chemtools/taxonomy/data/new_system.md`):
  - `calculable_features.json` (foundation): atomic SMARTS/count detectors
  - `calculable_features_properties.json`: heuristic/descriptor features
  - `calculable_features_derived.json`: boolean composition + derived_shortcuts

Callers should use `load_calculable_feature_spec()` to obtain the merged spec
in the legacy single-dict format expected by `chemtools.featurizers.calculable`.
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
    """Return the taxonomy data directory root for calculable feature assets."""
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


@lru_cache(maxsize=4)
def load_calculable_feature_spec(root: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load and merge the layered calculable feature specification.

    Args:
        root: Optional path to a taxonomy data directory containing the layered
              JSON files. Defaults to `chemtools/taxonomy/data`.

    Returns:
        A merged spec dict with keys: version, description, schema_notes,
        features, derived_shortcuts, changelog (when present).
    """
    data_root = _resolve_data_root(root)

    base_path = data_root / _BASE_FILE
    if not base_path.exists():
        raise FileNotFoundError(f"Missing calculable feature base spec: {base_path}")

    base = _safe_load_json(base_path)
    merged: Dict[str, Any] = dict(base)
    merged["features"] = _coerce_feature_list(base.get("features"), path=base_path)

    # Start with any derived_shortcuts in the base spec (expected to be absent in v6 layered mode).
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

    # Ensure legacy callers always see list types.
    if "changelog" in merged and merged["changelog"] is None:
        merged["changelog"] = []

    return merged


__all__ = ["load_calculable_feature_spec"]

