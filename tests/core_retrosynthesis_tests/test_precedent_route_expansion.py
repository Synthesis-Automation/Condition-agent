"""Two-step precedent-route chemical-space expansion regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

from core_retrosynthesis.precedent_route_expansion import (
    EXPANSION_LEVELS,
    PrecedentRouteExpansionInput,
    expand_two_step_precedent_route,
    load_precedent_route_expansion_definition,
    run_precedent_route_expansion_poc,
)


@pytest.fixture(scope="module")
def definition():
    return load_precedent_route_expansion_definition()


@pytest.fixture(scope="module")
def amine_result(definition):
    return expand_two_step_precedent_route(
        definition.routes[0], definition_id=definition.definition_id
    )


def test_default_panel_contains_three_distinct_two_step_routes(definition) -> None:
    assert definition.definition_id == (
        "two_step_precedent_route_expansion_poc.v1@1.0"
    )
    assert [route.route_id for route in definition.routes] == [
        "two_step_amine_acylation",
        "two_step_alcohol_acylation",
        "two_step_thiol_alkylation",
    ]


def test_expansion_replays_exact_route_and_is_strictly_nested(
    amine_result,
) -> None:
    assert amine_result.exact_replay_valid
    assert amine_result.warnings == ()
    assert [level.level for level in amine_result.levels] == list(
        EXPANSION_LEVELS
    )
    assert [level.product_count for level in amine_result.levels] == [1, 9, 20]
    assert [level.new_product_count for level in amine_result.levels] == [1, 8, 11]
    products = [set(level.product_smiles) for level in amine_result.levels]
    assert products[0] < products[1] < products[2]
    assert amine_result.target_smiles in products[0]
    assert "CC(=O)Nc1ccccc1" not in products[1]
    assert "CC(=O)Nc1ccccc1" in products[2]


def test_every_pathway_retains_forward_and_reverse_graph_evidence(
    amine_result,
) -> None:
    for pathway in amine_result.levels[-1].pathways:
        assert pathway.first_reverse_round_trip
        assert pathway.first_operator_edit_agreement
        assert pathway.second_reverse_round_trip
        assert pathway.second_operator_edit_agreement
        assert pathway.first_operator_id in amine_result.first_step_operator_ids
        assert pathway.second_operator_id in amine_result.second_step_operator_ids


def test_invalid_declared_input_is_counted_without_weakening_valid_space(
    definition,
) -> None:
    route = definition.routes[0]
    expanded = replace(
        route,
        first_step_inputs=(
            *route.first_step_inputs,
            PrecedentRouteExpansionInput("not-smiles", "R2"),
        ),
    )

    result = expand_two_step_precedent_route(
        expanded, definition_id=definition.definition_id
    )

    assert result.exact_replay_valid
    assert result.levels[-1].rejection_counts == {
        "first_step_no_supported_product": 1
    }


def test_panel_report_is_deterministic_and_serializable(
    definition,
    tmp_path,
) -> None:
    first_path = tmp_path / "first.json"
    second_path = tmp_path / "second.json"

    first = run_precedent_route_expansion_poc(
        definition, output_path=first_path
    )
    second = run_precedent_route_expansion_poc(
        definition, output_path=second_path
    )

    assert first == second
    assert first_path.read_bytes() == second_path.read_bytes()
    assert first["route_count"] == 3
    assert first["exact_replay_count"] == 3
    assert first["level_product_counts"] == {"R0": 3, "R1": 24, "R2": 55}
