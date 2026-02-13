"""Tests for v1-v2 benchmark utility."""

from poc_gpt52_reaction_v2.benchmark import run_benchmark


def test_benchmark_reports_consistent_summary_shape() -> None:
    report = run_benchmark()
    summary = report["summary"]

    assert summary["n_cases"] >= 5
    assert 0.0 <= summary["v1_mapped_exact_accuracy"] <= 1.0
    assert 0.0 <= summary["v2_exact_accuracy"] <= 1.0
    assert summary["v2_exact_accuracy"] >= summary["v1_mapped_exact_accuracy"]
    assert len(report["cases"]) == summary["n_cases"]
