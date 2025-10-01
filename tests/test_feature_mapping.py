from __future__ import annotations

from chemtools.agent.features.mapping import build_features, map_electrophile, map_nucleophile


def test_map_electrophile_detects_aryl_chloride() -> None:
    info = map_electrophile(["Clc1ccc(Cl)cc1"])
    assert info["class"] == "aryl chloride"


def test_map_nucleophile_detects_primary_aniline() -> None:
    info = map_nucleophile(["Nc1ccccc1"])
    assert info["class"] in {"primary aniline", "primary amine"}


def test_build_features_merges_context() -> None:
    features = build_features(["Clc1ccc(Cl)cc1", "Nc1ccccc1"], context={"mode": "batch", "temperature_C": 100})
    assert features["electrophile"]["class"].startswith("aryl")
    assert features["nucleophile"]["class"].startswith("primary")
    assert features["mode"] == "batch"
    assert features["temperature_C"] == 100
