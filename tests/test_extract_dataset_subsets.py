from __future__ import annotations

import csv
import json
from pathlib import Path

from scripts import extract_dataset_subsets as subsets


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reaction_id", "reaction_smiles", "reference"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_build_dataset_subsets_creates_per_file_outputs(tmp_path: Path) -> None:
    input_dir = tmp_path / "in"
    output_dir = tmp_path / "out"

    _write_csv(
        input_dir / "a.csv",
        [
            {"reaction_id": "a1", "reaction_smiles": "A>>B", "reference": "r1"},
            {"reaction_id": "a2", "reaction_smiles": "C>>D", "reference": "r1"},
            {"reaction_id": "a3", "reaction_smiles": "E>>F", "reference": "r2"},
            {"reaction_id": "a4", "reaction_smiles": "A>>B", "reference": "r3"},
            {"reaction_id": "a5", "reaction_smiles": "bad", "reference": "r4"},
        ],
    )
    _write_csv(
        input_dir / "b.csv",
        [
            {"reaction_id": "b1", "reaction_smiles": "G>>H", "reference": "r1"},
            {"reaction_id": "b2", "reaction_smiles": "I>>J", "reference": "r2"},
        ],
    )

    summary = subsets.build_dataset_subsets(
        input_dir=input_dir,
        output_dir=output_dir,
        per_file=2,
        seed=20260211,
        reaction_column="reaction_smiles",
        require_arrow=True,
    )

    assert summary["files_scanned"] == 2
    assert summary["total_sampled_rows"] == 4
    assert (output_dir / "a.csv").exists()
    assert (output_dir / "b.csv").exists()

    with (output_dir / "a.csv").open("r", encoding="utf-8", newline="") as handle:
        rows_a = list(csv.DictReader(handle))
    with (output_dir / "b.csv").open("r", encoding="utf-8", newline="") as handle:
        rows_b = list(csv.DictReader(handle))

    assert len(rows_a) == 2
    assert len(rows_b) == 2
    assert all(">>" in str(row["reaction_smiles"]) for row in rows_a + rows_b)

    # Ensure JSON-serializable summary contract for downstream tooling.
    json.dumps(summary)
