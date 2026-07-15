import csv
import json
from pathlib import Path

from reactive_taxonomy import (
    evaluate_molecular_features,
    load_molecular_feature_benchmark,
)


def test_phase1_benchmark_is_versioned_and_has_unique_cases() -> None:
    benchmark = load_molecular_feature_benchmark()
    assert benchmark["benchmark_id"] == "molecular_features.v1"
    assert len(benchmark["cases"]) >= 20
    case_ids = [case["case_id"] for case in benchmark["cases"]]
    assert len(case_ids) == len(set(case_ids))
    assert {case["case_id"] for case in benchmark["cases"]} >= {
        "secondary_dimethylamine",
        "primary_tert_butylamine",
        "hindered_electron_poor_aryl_bromide",
        "negative_ethane",
        "invalid_smiles",
    }
    partition_ids = [
        case_id
        for case_ids in benchmark["partitions"].values()
        for case_id in case_ids
    ]
    assert sorted(partition_ids) == sorted(case_ids)
    assert set(benchmark["partitions"]) == {"development", "validation", "test"}


def test_phase1_machine_gate_and_review_artifacts(tmp_path: Path) -> None:
    report = evaluate_molecular_features(tmp_path)
    assert report["machine_gate_passed"]
    assert report["human_gate_status"] == "pending_chemist_review"
    assert report["failed_case_ids"] == []
    metrics = report["metrics"]
    assert metrics["critical_case_pass_rate"] == 1.0
    assert metrics["functional_group_precision"] == 1.0
    assert metrics["reactive_site_recall"] == 1.0
    assert metrics["reactive_atom_accuracy"] == 1.0
    assert metrics["reactive_atom_check_count"] >= 20
    assert metrics["environment_check_accuracy"] == 1.0
    assert metrics["determinism_rate"] == 1.0
    assert metrics["equivalent_smiles_invariance_rate"] == 1.0
    assert metrics["taxonomy_error_count"] == 0
    assert set(report["partition_metrics"]) == {"development", "validation", "test"}
    assert all(
        partition["case_pass_rate"] == 1.0
        for partition in report["partition_metrics"].values()
    )
    expected_files = {
        "machine_report.json",
        "case_results.jsonl",
        "chemist_review.csv",
        "review_structures.html",
        "disagreements.csv",
    }
    assert expected_files <= {path.name for path in tmp_path.iterdir()}


def test_chemist_packet_is_blind_and_has_review_fields(tmp_path: Path) -> None:
    report = evaluate_molecular_features(tmp_path)
    with (tmp_path / "chemist_review.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == report["metrics"]["case_count"]
    assert "expected_group_counts" not in rows[0]
    assert "expected_site_counts" not in rows[0]
    for field in (
        "reviewer_id",
        "reactive_sites_correct",
        "functional_groups_correct",
        "steric_correct",
        "electronic_context_correct",
        "missing_features",
        "incorrect_features",
        "comments",
    ):
        assert field in rows[0]
        assert all(row[field] == "" for row in rows)
    review_html = (tmp_path / "review_structures.html").read_text(
        encoding="utf-8"
    )
    assert review_html.count("<section>") == len(rows)
    assert "benchmark expectations" in review_html
    assert "detected reactive-site atoms" in review_html


def test_phase1_evaluation_is_deterministic(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first_report = evaluate_molecular_features(first)
    second_report = evaluate_molecular_features(second)
    assert first_report == second_report
    assert (first / "case_results.jsonl").read_bytes() == (
        second / "case_results.jsonl"
    ).read_bytes()
    assert (first / "review_structures.html").read_bytes() == (
        second / "review_structures.html"
    ).read_bytes()


def test_case_results_retain_machine_failures_separately(tmp_path: Path) -> None:
    evaluate_molecular_features(tmp_path)
    cases = [
        json.loads(line)
        for line in (tmp_path / "case_results.jsonl").read_text(
            encoding="utf-8"
        ).splitlines()
    ]
    assert cases
    assert all(case["case_passed"] for case in cases)
    invalid = next(case for case in cases if case["case_id"] == "invalid_smiles")
    assert not invalid["valid"]
    assert invalid["error"] == "INVALID_SMILES"
