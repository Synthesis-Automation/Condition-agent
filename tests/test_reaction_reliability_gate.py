from __future__ import annotations

from scripts import reaction_reliability_gate as gate


def test_evaluate_gate_passes_without_baseline() -> None:
    report = {
        "unknown_reaction_type": {"rate": 0.2},
        "low_reaction_key_quality": {"rate": 0.01},
        "empty_reaction_key": {"rate": 0.0},
    }
    passed, errors, details = gate.evaluate_gate(
        report,
        max_unknown_rate=0.3,
        max_low_quality_rate=0.02,
        max_empty_key_rate=0.01,
    )
    assert passed is True
    assert not errors
    assert details["current"]["unknown_reaction_type"] == 0.2


def test_evaluate_gate_fails_on_threshold_and_delta() -> None:
    current = {
        "unknown_reaction_type": {"rate": 0.5},
        "low_reaction_key_quality": {"rate": 0.03},
        "empty_reaction_key": {"rate": 0.02},
    }
    baseline = {
        "unknown_reaction_type": {"rate": 0.3},
        "low_reaction_key_quality": {"rate": 0.01},
        "empty_reaction_key": {"rate": 0.0},
    }
    passed, errors, details = gate.evaluate_gate(
        current,
        max_unknown_rate=0.4,
        max_low_quality_rate=0.02,
        max_empty_key_rate=0.01,
        baseline=baseline,
        max_unknown_delta=0.05,
        max_low_quality_delta=0.01,
        max_empty_key_delta=0.005,
    )
    assert passed is False
    assert any(err.startswith("threshold_failed:unknown_reaction_type") for err in errors)
    assert any(err.startswith("delta_failed:unknown_reaction_type") for err in errors)
    assert details["deltas"]["unknown_reaction_type"] == 0.2
