"""
Lightweight tests for the constraint parser utilities.

These tests run without touching the LangChain runtime and validate that
the deterministic parsing helpers behave as expected.
"""

from lang_chain.constraint_parser import (
    build_constraint_spec,
    filter_cores_by_constraints,
)


def test_build_constraint_spec_from_text():
    spec = build_constraint_spec(text="Pd-free, prefer copper, search all families")
    assert "PD" in spec.exclude_metals
    assert "CU" in spec.prefer_metals
    assert spec.search_all_families is True


def test_build_constraint_spec_from_lists():
    spec = build_constraint_spec(
        allow_metals=["Ni", "pd"],
        exclude_metals=["Cu"],
        prefer_metals=["Ni"],
    )
    assert spec.allow_metals == {"NI", "PD"}
    assert spec.exclude_metals == {"CU"}
    assert spec.prefer_metals == {"NI"}


def test_filter_cores_by_constraints():
    cores = ["Pd/XPhos", "Cu/phen", "Ni/Bpy"]
    filtered = filter_cores_by_constraints(
        cores,
        allow_metals={"CU", "NI"},
        exclude_metals={"CU"},
        prefer_metals={"NI"},
    )
    assert filtered[0] == "Ni/Bpy"
    assert "Pd/XPhos" not in filtered
