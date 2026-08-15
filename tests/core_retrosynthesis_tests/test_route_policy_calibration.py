"""Whole-route route-policy calibration regressions."""

from __future__ import annotations

import pytest

from core_retrosynthesis.route_action_evaluation import RouteActionEvaluation
from core_retrosynthesis.route_policy_calibration import (
    RoutePolicyCalibrationConfig,
    RoutePolicyScaleOutcome,
    build_route_policy_calibration_targets,
    select_route_policy_scale,
)


def _outcome(
    scale: float,
    *,
    solved: int = 0,
    strategies: int = 0,
    terminal_fraction: float = 0.0,
    expansions: int = 10,
) -> RoutePolicyScaleOutcome:
    return RoutePolicyScaleOutcome(
        residual_scale=scale,
        target_count=2,
        solved_target_count=solved,
        observed_strategy_recovered_count=strategies,
        mean_best_terminal_fraction=terminal_fraction,
        first_solution_efficiency=0,
        expanded_state_count=expansions,
        policy_reordered_expansion_count=0,
        targets=(),
    )


def _evaluation(*, split: str = "validation") -> RouteActionEvaluation:
    return RouteActionEvaluation(
        evaluation_id="evaluation",
        tree_id="tree",
        source_route_id="route",
        patent_id="US1A1",
        split=split,
        target_smiles="CCCC",
        reaction_count=0,
        maximum_depth=3,
        search_config_id="search",
        steps=(),
    )


def test_scale_selection_prefers_zero_when_route_metrics_tie() -> None:
    selected = select_route_policy_scale(
        (_outcome(0.0), _outcome(0.25), _outcome(1.0))
    )

    assert selected.residual_scale == 0.0


def test_scale_selection_uses_route_completion_before_policy_strength() -> None:
    selected = select_route_policy_scale(
        (
            _outcome(0.0, solved=1, strategies=2),
            _outcome(0.25, solved=1, strategies=3),
            _outcome(1.0, solved=2, strategies=1),
        )
    )

    assert selected.residual_scale == 1.0


def test_calibration_targets_must_come_from_requested_split() -> None:
    targets = build_route_policy_calibration_targets(
        (_evaluation(),), ("route",)
    )

    assert targets[0].route_id == "route"
    assert targets[0].target_smiles == "CCCC"

    with pytest.raises(ValueError, match="not in validation split"):
        build_route_policy_calibration_targets(
            (_evaluation(split="test"),), ("route",)
        )


def test_calibration_config_is_versioned_and_validates_budget() -> None:
    config = RoutePolicyCalibrationConfig()

    assert config.config_id.startswith("RPCC1:")
    assert config.planner_kwargs()["max_expansions"] == 15
    assert "minimum_validation_targets" not in config.planner_kwargs()

    with pytest.raises(ValueError, match="cover per-step top-k"):
        RoutePolicyCalibrationConfig(
            per_step_top_k=5,
            max_candidates_to_validate=4,
        )
