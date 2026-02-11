from __future__ import annotations

import csv
from pathlib import Path

from scripts import extract_reference_testset as sampler


def _write_csv(path: Path) -> None:
    rows = [
        {"reaction_id": "r1", "reaction_smiles": "A>>B", "reference": "Ref-A"},
        {"reaction_id": "r2", "reaction_smiles": "C>>D", "reference": "Ref-A"},
        {"reaction_id": "r3", "reaction_smiles": "E>>F", "reference": "Ref-A"},
        {"reaction_id": "r4", "reaction_smiles": "G>>H", "reference": "Ref-B"},
        {"reaction_id": "r5", "reaction_smiles": "I>>J", "reference": "Ref-B"},
        {"reaction_id": "r6", "reaction_smiles": "bad", "reference": "Ref-C"},
        {"reaction_id": "r7", "reaction_smiles": "K>>L", "reference": ""},
        {"reaction_id": "r8", "reaction_smiles": "A>>B", "reference": "Ref-A"},
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reaction_id", "reaction_smiles", "reference"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_build_reference_samples_picks_two_per_reference(tmp_path: Path) -> None:
    input_csv = tmp_path / "input.csv"
    _write_csv(input_csv)

    rows, summary, fieldnames = sampler.build_reference_samples(
        input_csv=input_csv,
        reference_column="reference",
        reaction_column="reaction_smiles",
        per_reference=2,
        seed=20260211,
        require_arrow=True,
    )

    assert list(fieldnames) == ["reaction_id", "reaction_smiles", "reference"]
    assert summary["unique_references"] == 2
    assert summary["sampled_rows_total"] == 4
    assert summary["references_with_lt_target"] == 0

    by_ref = {}
    for row in rows:
        ref = str(row["reference"])
        by_ref.setdefault(ref, []).append(row)

    assert set(by_ref.keys()) == {"Ref-A", "Ref-B"}
    assert len(by_ref["Ref-A"]) == 2
    assert len(by_ref["Ref-B"]) == 2
    assert all(">>" in str(r["reaction_smiles"]) for r in rows)
    assert all(str(r["sample_rank_in_reference"]) in {"1", "2"} or isinstance(r["sample_rank_in_reference"], int) for r in rows)
