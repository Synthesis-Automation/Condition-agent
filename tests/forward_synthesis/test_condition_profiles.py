"""Condition-profile definition and validation regressions."""

from __future__ import annotations

import pytest

from forward_synthesis import (
    assess_condition_profile,
    condition_profile_catalog,
    normalize_condition_profile,
)


def test_condition_profile_catalog_is_versioned_and_complete() -> None:
    catalog = condition_profile_catalog()

    assert catalog["definition_id"] == "forward_condition_profiles.v1"
    assert {item["id"] for item in catalog["strategies"]} >= {
        "unspecified",
        "transition_metal_catalysis",
        "metal_free_polar",
        "radical",
        "photochemical",
        "thermal",
    }


def test_catalyst_family_requires_transition_metal_strategy() -> None:
    with pytest.raises(ValueError, match="requires transition_metal"):
        normalize_condition_profile(
            {"strategy": "thermal", "catalyst_family": "palladium"}
        )


def test_redox_profiles_use_net_edit_evidence_without_hard_filtering() -> None:
    reduction = assess_condition_profile(
        (
            "broken:Br-C:SINGLE>NONE",
            "hydrogen_change:C-H:NONE>SINGLE",
        ),
        {"redox_mode": "reductive"},
    )
    oxidation = assess_condition_profile(
        ("hydrogen_change:C-H:SINGLE>NONE",),
        {"redox_mode": "oxidative"},
    )

    assert reduction.score_adjustment == 0.08
    assert reduction.matched_rules == ("reductive_pattern_match",)
    assert oxidation.score_adjustment == 0.08
    assert oxidation.matched_rules == ("oxidative_pattern_match",)
