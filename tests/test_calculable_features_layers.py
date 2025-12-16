"""
Regression tests for the layered calculable feature spec.

The spec is split into:
  - calculable_features.json (foundation)
  - calculable_features_properties.json
  - calculable_features_derived.json

This test ensures the runtime loader merges layers correctly and enforces the
intended separation between atomic SMARTS, heuristic properties, and derived logic.
"""

from __future__ import annotations

import json
from pathlib import Path

from chemtools.taxonomy.calculable_spec import load_calculable_feature_spec


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_layer_files_exist_and_are_disjoint():
    data_root = Path(__file__).resolve().parent.parent / "chemtools" / "taxonomy" / "data"
    base = _load(data_root / "calculable_features.json")
    props = _load(data_root / "calculable_features_properties.json")
    derived = _load(data_root / "calculable_features_derived.json")

    base_features = base.get("features") or []
    props_features = props.get("features") or []
    derived_features = derived.get("features") or []

    base_tokens = {f.get("token") for f in base_features if isinstance(f, dict) and f.get("token")}
    props_tokens = {f.get("token") for f in props_features if isinstance(f, dict) and f.get("token")}
    derived_tokens = {f.get("token") for f in derived_features if isinstance(f, dict) and f.get("token")}

    assert base_tokens
    assert props_tokens
    assert derived_tokens

    assert not (base_tokens & props_tokens)
    assert not (base_tokens & derived_tokens)
    assert not (props_tokens & derived_tokens)


def test_foundation_is_atomic_and_nonheuristic():
    data_root = Path(__file__).resolve().parent.parent / "chemtools" / "taxonomy" / "data"
    base = _load(data_root / "calculable_features.json")

    assert "derived_shortcuts" not in base
    features = base.get("features") or []

    for feat in features:
        if not isinstance(feat, dict):
            continue
        expr = feat.get("derive") or feat.get("derived")
        assert not (isinstance(expr, str) and expr.strip()), feat.get("token")
        detect = feat.get("detect")
        assert not (isinstance(detect, dict) and "heuristic" in detect), feat.get("token")


def test_loader_merges_layers():
    data_root = Path(__file__).resolve().parent.parent / "chemtools" / "taxonomy" / "data"
    base = _load(data_root / "calculable_features.json")
    props = _load(data_root / "calculable_features_properties.json")
    derived = _load(data_root / "calculable_features_derived.json")

    merged = load_calculable_feature_spec()
    assert merged.get("version") == base.get("version")

    expected_feature_count = len(base.get("features") or []) + len(props.get("features") or []) + len(derived.get("features") or [])
    assert len(merged.get("features") or []) == expected_feature_count
    assert len(merged.get("derived_shortcuts") or []) == len(derived.get("derived_shortcuts") or [])

