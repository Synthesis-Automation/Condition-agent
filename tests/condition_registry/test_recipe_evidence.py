"""Agreement, explicit absence and missing operating data in recipe queries."""

import pytest

from condition_registry.recipe_evidence import compare_evidence_recipes, normalize_evidence_recipe


def component(name, role, **quantity):
    return {"raw_identifier": name, "source_role_hint": role, **quantity}


BASE = component("cesium carbonate", "base", amount=1.5, amount_unit="equiv")
SOLVENT = component("DMSO", "solvent")


def compare(request, source):
    return compare_evidence_recipes(normalize_evidence_recipe(request), normalize_evidence_recipe(source))


def test_registry_aliases_quantities_and_stage_operating_values():
    query = {"components": [BASE, SOLVENT], "temperature_c": 25, "time_h": 1,
             "illumination": "on", "atmosphere": "nitrogen"}
    source = {"components": [component("67-68-5", "solvent"), BASE],
              "stages": [{"stage_index": 1, "time_h": 1, "temperature_c": 25,
                          "atmosphere": "N2", "provenance": {"illumination": "white LEDs"}}]}
    comparison = compare(query, source)
    assert comparison.status == "matches_stated_fields" and comparison.source_only
    assert any("1.5" in v for v in comparison.agreements)
    assert compare({**query, "illumination": "blue LEDs"}, source).status == "different"
    assert compare({**query, "components": [{**BASE, "amount": 0.5}]}, source).status == "different"
    assert compare({**query, "components": [{**BASE, "amount_unit": "mg"}]}, source).status == "partially_reported"


def test_omissions_require_source_evidence_and_missing_light_is_unknown():
    query = {"components": [SOLVENT], "declared_absences": ["base"], "illumination": False}
    assert compare(query, {"components": [SOLVENT]}).status == "partially_reported"
    comparison = compare(query, {"components": [BASE, SOLVENT], "illumination": True})
    assert comparison.status == "different" and len(comparison.differences) == 2
    assert compare(query, {"components": [SOLVENT], "declared_absences": ["base", "light"]}).status == "matches_stated_fields"
    assert compare({"declared_absences": ["base"]}, {"time_h": 1}).status == "unknown"


def test_multistage_times_and_different_stages_are_not_inferred():
    source = {"stages": [{"stage_index": 1, "time_h": 1, "temperature_c": 25},
                         {"stage_index": 2, "time_h": 1, "temperature_c": 50}]}
    assert compare({"time_h": 1}, source).status == "unknown"
    assert compare({"temperature_c": 25}, source).status == "unknown"
    assert compare(source, source).status == "matches_stated_fields"
    reversed_rows = {"stages": list(reversed(source["stages"]))}
    assert compare(source, reversed_rows) == compare(source, source)
    unknown = {"components": [component("yet-unidentified-condition", "base")]}
    assert compare(unknown, unknown).status == "unknown"


def test_known_ingredient_and_amount_survive_unassigned_solvent_role():
    carrier = component("CCc1cc(CC)cc(CC)c1", "solvent", identifier_type="smiles", amount=1, amount_unit="uL")
    query = {"components": [BASE, carrier]}
    comparable = compare(query, query)
    assert comparable.status == "partially_reported"
    assert any(v.startswith("ingredient/cas:102-25-0") for v in comparable.agreements)
    assert any(v.startswith("amount/cas:102-25-0") for v in comparable.agreements)
    assert any(v.startswith("reported_role_hint/cas:102-25-0") for v in comparable.agreements)
    assert any("unassigned" in v for v in comparable.unknowns)
    alternate = compare(query, {"components": [component("sodium hydroxide", "base"), carrier]})
    assert alternate.status == "different"
    assert any("bases/cas:" in v for v in alternate.differences)


def test_unknown_illumination_and_conflicting_dark_stages_cannot_support_a_recipe():
    assert compare({"illumination": "on"}, {"illumination": "unreported"}).status == "unknown"
    assert compare({"illumination": True}, {"illumination": "unrecognized entry"}).status == "unknown"
    dark = normalize_evidence_recipe({"declared_absences": ["light"]})
    conflicting = {"declared_absences": ["light"], "stages": [
        {"stage_index": 1, "provenance": {"illumination": "off"}},
        {"stage_index": 2, "provenance": {"illumination": "white LEDs"}}]}
    assert compare_evidence_recipes(dark, conflicting).status == "unknown"
    with pytest.raises(ValueError, match="Contradictory"):
        normalize_evidence_recipe(conflicting)


def test_zero_recorded_component_amount_does_not_confirm_presence_or_omission():
    query = normalize_evidence_recipe({"components": [component("cesium carbonate", "base")]})
    source = normalize_evidence_recipe({"components": [BASE]})
    source["bases"][0]["amount"] = 0
    assert compare_evidence_recipes(query, source).status == "unknown"
    assert compare_evidence_recipes(normalize_evidence_recipe({"declared_absences": ["base"]}), source).status == "unknown"
    with pytest.raises(ValueError, match="declared_absences"):
        normalize_evidence_recipe({"components": [{**BASE, "amount": 0}]})


@pytest.mark.parametrize("bad", [
    {}, {"temperature_c": float("nan")}, {"typo": 2},
    {"components": [BASE], "declared_absences": ["base"]},
    {"illumination": "white LEDs", "declared_absences": ["light"]},
    {"components": [{**BASE, "identifier_type": []}]},
    {"components": [{**BASE, "amount_unit": None}]},
    {"components": [{**BASE, "amount": -1}]},
    {"recipe_id": "r", "bases": [{"substance_id": "x", "amount": float("inf"), "amount_unit": "g"}]},
    {"temperature_c": 25, "stages": [{"stage_index": 1, "temperature_c": 45}]},
])
def test_invalid_or_contradictory_queries_fail(bad):
    with pytest.raises(ValueError):
        normalize_evidence_recipe(bad)
