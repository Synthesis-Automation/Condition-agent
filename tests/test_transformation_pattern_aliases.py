from __future__ import annotations

from chemtools.featurizers.formatters.aggregation import (
    _resolve_transformation_pattern_info,
    load_transformation_patterns,
)


def test_click_canonical_taxonomy_id_resolves_to_transformation_pattern() -> None:
    payload = load_transformation_patterns()
    info = _resolve_transformation_pattern_info(
        "Click_azide_alkyne_cycloaddition",
        payload,
    )
    assert info is not None
    assert info.get("pattern") == "cycloaddition"


def test_wacker_canonical_taxonomy_id_resolves_split_transformation_patterns() -> None:
    payload = load_transformation_patterns()
    info = _resolve_transformation_pattern_info("Wacker_oxidation", payload)
    assert info is not None
    assert info.get("pattern") == "addition"
    # Split variants should merge into an effective rule payload.
    assert "leaving_groups" in info
