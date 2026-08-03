"""Batch-audit regressions for display-minimized reaction graphics."""

from __future__ import annotations

import csv

import pytest

from benchmarks.reaction_display_batch_audit import build_batch_audit


def _write_source(path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "canonical_reaction_smiles",
                "reaction_label",
                "original_reaction_type",
                "reaction_core_status",
                "reaction_core_quality_status",
            ),
        )
        writer.writeheader()
        writer.writerows(
            (
                {
                    "canonical_reaction_smiles": "CCBr.N>>CCN",
                    "reaction_label": "test coupling",
                },
                {
                    "canonical_reaction_smiles": "CC.O>>N#N",
                    "reaction_label": "unsupported accounting",
                    "reaction_core_status": "unavailable",
                },
                {
                    "canonical_reaction_smiles": "CCBr.N>>CCN",
                    "reaction_label": "duplicate coupling",
                },
            )
        )


def test_batch_audit_writes_one_graph_per_successful_row_and_failures(
    tmp_path,
) -> None:
    source = tmp_path / "reaction_review.csv"
    output = tmp_path / "audit"
    _write_source(source)

    summary = build_batch_audit(
        source_path=source,
        output_dir=output,
        use_rxnmapper=False,
    )

    assert summary["input_rows"] == 3
    assert summary["unique_reactions"] == 2
    assert summary["successful_rows"] == 2
    assert summary["unsuccessful_rows"] == 1
    assert summary["failure_by_stage"] == {"reaction_core": 1}
    graphs = sorted((output / "graphs").glob("*.svg"))
    assert len(graphs) == 2
    assert all(path.read_bytes().startswith(b"<?xml") for path in graphs)
    with (output / "unsuccessful_rows.csv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        failures = list(csv.DictReader(handle))
    assert len(failures) == 1
    assert failures[0]["source_data_row"] == "2"
    assert failures[0]["source_csv_line"] == "3"
    assert failures[0]["failure_stage"] == "reaction_core"
    assert failures[0]["canonical_reaction_smiles"] == "CC.O>>N#N"


def test_batch_audit_refuses_nonempty_output_directory(tmp_path) -> None:
    source = tmp_path / "reaction_review.csv"
    output = tmp_path / "audit"
    _write_source(source)
    output.mkdir()
    (output / "existing.txt").write_text("preserve", encoding="utf-8")

    with pytest.raises(FileExistsError, match="output directory is not empty"):
        build_batch_audit(
            source_path=source,
            output_dir=output,
            use_rxnmapper=False,
        )

    assert (output / "existing.txt").read_text(encoding="utf-8") == "preserve"
