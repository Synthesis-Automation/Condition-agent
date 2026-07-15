import csv
import json
from pathlib import Path

from reactive_taxonomy import evaluate_reaction_edits, load_reaction_edit_benchmark


def test_phase2_benchmark_is_versioned_and_partitioned() -> None:
    benchmark = load_reaction_edit_benchmark()
    assert benchmark["benchmark_id"] == "reaction_edits.v1"
    assert len(benchmark["cases"]) >= 12
    case_ids = [case["case_id"] for case in benchmark["cases"]]
    partition_ids = [
        case_id
        for ids_in_partition in benchmark["partitions"].values()
        for case_id in ids_in_partition
    ]
    assert len(case_ids) == len(set(case_ids))
    assert sorted(partition_ids) == sorted(case_ids)
    assert set(benchmark["partitions"]) == {"development", "validation", "test"}


def test_phase2_machine_gate_and_artifacts(tmp_path: Path) -> None:
    report = evaluate_reaction_edits(tmp_path)

    assert report["machine_gate_passed"]
    assert report["human_gate_status"] == "pending_chemist_review"
    assert report["failed_case_ids"] == []
    assert all(report["threshold_results"].values())
    assert report["metrics"]["edit_precision"] == 1.0
    assert report["metrics"]["edit_recall"] == 1.0
    assert report["metrics"]["mapped_unmapped_parity_rate"] == 1.0
    assert report["metrics"]["conflict_retention_rate"] == 1.0
    assert {
        "machine_report.json",
        "case_results.jsonl",
        "chemist_review.csv",
        "review_structures.html",
        "disagreements.csv",
    } <= {path.name for path in tmp_path.iterdir()}


def test_phase2_review_packet_is_blind_and_readable(tmp_path: Path) -> None:
    report = evaluate_reaction_edits(tmp_path)
    with (tmp_path / "chemist_review.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == report["metrics"]["case_count"]
    assert "expected_edits" not in rows[0]
    assert "expected_evidence" not in rows[0]
    for field in (
        "reviewer_id",
        "reaction_center_correct",
        "bond_edits_correct",
        "hydrogen_changes_correct",
        "evidence_assessment_correct",
        "missing_edits",
        "incorrect_edits",
        "comments",
    ):
        assert field in rows[0]
        assert all(row[field] == "" for row in rows)
    hydrogenation = next(
        row for row in rows if row["case_id"] == "mapped_alkene_hydrogenation"
    )
    assert "change C:1–C:2: double→single" in hydrogenation["edit_descriptions"]
    assert "H gain at C:1" in hydrogenation["edit_descriptions"]
    review_html = (tmp_path / "review_structures.html").read_text(encoding="utf-8")
    assert review_html.count("<section>") == len(rows)
    assert "benchmark expectations" in review_html
    assert "reactant edit atoms" in review_html


def test_phase2_evaluation_is_deterministic(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first_report = evaluate_reaction_edits(first)
    second_report = evaluate_reaction_edits(second)

    assert first_report == second_report
    assert (first / "case_results.jsonl").read_bytes() == (
        second / "case_results.jsonl"
    ).read_bytes()
    assert (first / "review_structures.html").read_bytes() == (
        second / "review_structures.html"
    ).read_bytes()


def test_phase2_case_results_retain_conflicts_and_invalid_valence(
    tmp_path: Path,
) -> None:
    evaluate_reaction_edits(tmp_path)
    cases = [
        json.loads(line)
        for line in (tmp_path / "case_results.jsonl").read_text(
            encoding="utf-8"
        ).splitlines()
    ]

    conflict = next(
        case for case in cases if case["case_id"] == "mapping_reconstruction_conflict"
    )
    assert conflict["observed_evidence"] == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in conflict["warnings"]
    invalid = next(case for case in cases if case["case_id"] == "invalid_valence_product")
    assert not invalid["valid"]
    assert invalid["error"] == "INVALID_COMPONENT"
