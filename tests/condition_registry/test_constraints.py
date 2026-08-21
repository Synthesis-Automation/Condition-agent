"""Condition registry constraint normalization tests."""

from __future__ import annotations

from condition_registry import (
    ConditionConstraintSet,
    condition_constraint_conflicts,
    normalize_condition_constraint,
)


def test_substance_constraint_resolves_alias_to_stable_identity() -> None:
    resolution = normalize_condition_constraint(
        "excluded_substance",
        "Pd(PPh3)4",
        provenance="confirmed_user",
    )

    assert resolution.status == "resolved"
    assert resolution.constraint is not None
    assert resolution.constraint.normalized_value == "cas:14221-01-3"
    assert resolution.constraint.source_value == "Pd(PPh3)4"


def test_unknown_substance_and_atmosphere_are_not_guessed() -> None:
    substance = normalize_condition_constraint(
        "excluded_substance",
        "definitely-not-a-known-reagent",
        provenance="explicit_user",
    )
    atmosphere = normalize_condition_constraint(
        "required_atmosphere",
        "special inert mix",
        provenance="explicit_user",
    )

    assert substance.status == "unresolved"
    assert atmosphere.status == "unresolved"


def test_constraint_conflicts_use_canonical_recipe_fields() -> None:
    excluded = normalize_condition_constraint(
        "excluded_substance",
        "Pd(PPh3)4",
        provenance="confirmed_user",
    ).constraint
    temperature = normalize_condition_constraint(
        "maximum_temperature_c",
        60,
        provenance="explicit_user",
    ).constraint
    assert excluded is not None and temperature is not None
    recipe = {
        "catalysts": [
            {
                "substance_id": "cas:14221-01-3",
                "primary_role": "metal_catalyst",
            }
        ],
        "temperature_c": 80.0,
    }

    conflicts = condition_constraint_conflicts(
        recipe,
        ConditionConstraintSet((excluded, temperature)),
    )

    assert any(value.endswith(":excluded_substance") for value in conflicts)
    assert any(value.endswith(":temperature_exceeded") for value in conflicts)


def test_registry_rejects_non_user_hard_constraint_provenance() -> None:
    try:
        normalize_condition_constraint(
            "excluded_role",
            "metal_catalyst",
            provenance="model_proposed",  # type: ignore[arg-type]
        )
    except ValueError as exc:
        assert "user authority" in str(exc)
    else:
        raise AssertionError("model proposals cannot create hard constraints")
