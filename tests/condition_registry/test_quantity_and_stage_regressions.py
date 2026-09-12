"""Counterexamples for quantity identity and stage-wide user constraints."""

from dataclasses import replace

import pytest

from condition_registry import (
    ConditionConstraintSet,
    normalize_condition_constraint,
    condition_constraint_conflicts,
)
from condition_registry.models import ConditionComponentInput
from condition_registry.recipes import build_resolved_recipe_from_inputs


def _constraint(kind, value):
    resolved = normalize_condition_constraint(kind, value, provenance="explicit_user")
    assert resolved.constraint is not None
    return ConditionConstraintSet((resolved.constraint,))


@pytest.mark.parametrize(
    "stages",
    [
        [{"temperature_c": 120}],
        [{"temperature_c": None}],
        [{"temperature_c": 25}, {"temperature_c": 120}],
    ],
)
def test_maximum_temperature_checks_every_stage(stages):
    assert condition_constraint_conflicts(
        {"temperature_c": 25, "stages": stages},
        _constraint("maximum_temperature_c", 50),
    )


def test_stage_temperatures_can_supply_missing_summary():
    assert not condition_constraint_conflicts(
        {"stages": [{"temperature_c": 25}, {"temperature_c": 45}]},
        _constraint("maximum_temperature_c", 50),
    )


@pytest.mark.parametrize("atmosphere", ["air", "nitrogen-free", "argon/nitrogen", None])
def test_atmosphere_checks_stages_and_exact_gas_tokens(atmosphere):
    assert condition_constraint_conflicts(
        {"atmosphere": "N2", "stages": [{"atmosphere": atmosphere}]},
        _constraint("required_atmosphere", "nitrogen"),
    )


def test_non_finite_temperature_constraint_is_invalid():
    assert (
        normalize_condition_constraint(
            "maximum_temperature_c", "nan", provenance="explicit_user"
        ).status
        == "invalid"
    )


def _input(amount, unit="mL", **provenance):
    return ConditionComponentInput(
        raw_identifier="64-17-5",
        source_field="solvent_cas",
        amount=amount,
        amount_unit=unit,
        provenance=provenance,
    )


def test_conflicting_quantity_identity_is_order_invariant():
    values = [_input(1), _input(2)]
    recipes = [
        build_resolved_recipe_from_inputs(order)
        for order in (values, list(reversed(values)))
    ]
    assert recipes[0].to_dict() == recipes[1].to_dict()
    component = recipes[0].solvents[0]
    assert component.amount is None and component.quantity_status == "conflicting"
    assert len(component.quantity_observations) == 2
    assert "CONFLICTING_QUANTITY_OBSERVATIONS" in recipes[0].warnings


def test_equivalent_units_merge_without_summing_duplicate_reports():
    recipe = build_resolved_recipe_from_inputs([_input(1, "L"), _input(1000)])
    assert recipe.solvents[0].amount == 1000
    assert recipe.solvents[0].amount_unit == "mL"
    assert recipe.solvents[0].quantity_status == "reported"


def test_only_explicit_distinct_same_stage_additions_are_summed():
    values = [
        _input(
            amount,
            addition_id=str(number),
            stage_index=0,
            quantity_relationship="separate_addition",
        )
        for number, amount in enumerate((1, 2))
    ]
    assert build_resolved_recipe_from_inputs(values).solvents[0].amount == 3
    values[1] = replace(
        values[1], provenance={**values[1].provenance, "stage_index": 1}
    )
    assert build_resolved_recipe_from_inputs(values).solvents[0].amount is None
