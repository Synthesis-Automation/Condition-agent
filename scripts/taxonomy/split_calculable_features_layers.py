"""
Split `calculable_features.json` into 3 maintainable layers.

This implements the design philosophy in `chemtools/taxonomy/data/new_system.md`:
keep SMARTS as atomic structural sensors, and keep heavier logic in separate
layers.

Layers produced (in `chemtools/taxonomy/data/`):
  - `calculable_features.json` (foundation): SMARTS/count-based atomic features
  - `calculable_features_properties.json`: heuristic/descriptor-based features
  - `calculable_features_derived.json`: boolean derived rules + derived_shortcuts

Runtime code should load and merge these layers via `chemtools.taxonomy.calculable_spec`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple


def _classify_feature(feature: Dict[str, Any]) -> str:
    """
    Classify a feature dict into a layer.

    Rules:
      - Derived layer: has `derive`/`derived` string.
      - Properties layer: detect.heuristic present (and not derived).
      - Foundation layer: everything else (SMARTS/count atomic features).
    """
    derive_expr = feature.get("derive") or feature.get("derived")
    if isinstance(derive_expr, str) and derive_expr.strip():
        return "derived"

    detect = feature.get("detect")
    if isinstance(detect, dict) and "heuristic" in detect:
        return "properties"

    return "foundation"


def split_spec(payload: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Return (foundation_spec, properties_spec, derived_spec)."""
    features = payload.get("features", [])
    if not isinstance(features, list):
        raise ValueError("Expected top-level 'features' to be a list.")

    foundation_features: List[Dict[str, Any]] = []
    properties_features: List[Dict[str, Any]] = []
    derived_features: List[Dict[str, Any]] = []

    for entry in features:
        if not isinstance(entry, dict):
            continue
        layer = _classify_feature(entry)
        if layer == "derived":
            derived_features.append(entry)
        elif layer == "properties":
            properties_features.append(entry)
        else:
            foundation_features.append(entry)

    derived_shortcuts = payload.get("derived_shortcuts", [])
    if derived_shortcuts is None:
        derived_shortcuts = []
    if not isinstance(derived_shortcuts, list):
        raise ValueError("Expected top-level 'derived_shortcuts' to be a list.")

    # Foundation: keep metadata, drop derived_shortcuts.
    foundation_spec: Dict[str, Any] = {
        "version": payload.get("version"),
        "description": payload.get("description"),
        "schema_notes": payload.get("schema_notes"),
        "features": foundation_features,
        "changelog": payload.get("changelog", []),
    }

    version = str(payload.get("version") or "")
    properties_spec: Dict[str, Any] = {
        "version": version,
        "layer": "properties",
        "features": properties_features,
    }

    derived_spec: Dict[str, Any] = {
        "version": version,
        "layer": "derived",
        "features": derived_features,
        "derived_shortcuts": derived_shortcuts,
    }

    return foundation_spec, properties_spec, derived_spec


def _assert_disjoint_tokens(specs: List[Dict[str, Any]]) -> None:
    tokens_by_layer: Dict[str, set[str]] = {}
    for spec in specs:
        layer = str(spec.get("layer") or "foundation")
        tokens: set[str] = set()
        for feature in spec.get("features", []):
            if not isinstance(feature, dict):
                continue
            token = feature.get("token")
            if isinstance(token, str) and token:
                tokens.add(token)
        tokens_by_layer[layer] = tokens

    layers = list(tokens_by_layer.keys())
    for i, a in enumerate(layers):
        for b in layers[i + 1 :]:
            overlap = tokens_by_layer[a] & tokens_by_layer[b]
            if overlap:
                raise ValueError(f"Token overlap between {a} and {b}: {sorted(overlap)[:10]}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help=(
            "Optional path to a merged calculable feature spec. "
            "When omitted, loads the merged spec from the existing layered files via "
            "`chemtools.taxonomy.calculable_spec.load_calculable_feature_spec`."
        ),
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Optional taxonomy data directory root (defaults to parent of --input).",
    )
    args = parser.parse_args()

    root = args.root.resolve() if args.root is not None else None
    if root is None:
        if args.input is None:
            root = Path("chemtools/taxonomy/data").resolve()
        else:
            root = args.input.resolve().parent

    base_path = root / "calculable_features.json"
    properties_path = root / "calculable_features_properties.json"
    derived_path = root / "calculable_features_derived.json"

    if args.input is not None:
        payload = json.loads(args.input.resolve().read_text(encoding="utf-8"))
    else:
        from chemtools.archive.taxonomy.calculable_spec import load_calculable_feature_spec

        payload = load_calculable_feature_spec(root)
    foundation_spec, properties_spec, derived_spec = split_spec(payload)

    _assert_disjoint_tokens([foundation_spec, properties_spec, derived_spec])

    base_path.write_text(json.dumps(foundation_spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    properties_path.write_text(
        json.dumps(properties_spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    derived_path.write_text(json.dumps(derived_spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print("Wrote layered calculable feature specs:")
    print(f"- foundation: {base_path}")
    print(f"- properties: {properties_path}")
    print(f"- derived:    {derived_path}")
    print(f"Features: foundation={len(foundation_spec['features'])}, properties={len(properties_spec['features'])}, derived={len(derived_spec['features'])}")
    print(f"Derived shortcuts: {len(derived_spec.get('derived_shortcuts', []))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
