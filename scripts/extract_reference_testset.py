#!/usr/bin/env python3
"""
Build a deterministic test dataset by sampling N reactions per reference.

Example:
  python scripts/extract_reference_testset.py \
    --input data-processor/reaction_dataset/C_N_Coupling.csv \
    --output results/C_N_Coupling.reference2_test.csv \
    --summary-output results/C_N_Coupling.reference2_test.summary.json \
    --per-reference 2
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, List, Sequence, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


META_COLUMNS = [
    "source_file",
    "original_row_number",
    "sample_rank_in_reference",
    "rows_in_reference",
]


def _is_valid_reaction_smiles(value: Any, require_arrow: bool) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    if require_arrow and ">>" not in text:
        return False
    return True


def _clean_text(value: Any) -> str:
    return str(value or "").strip()


def build_reference_samples(
    *,
    input_csv: Path,
    reference_column: str,
    reaction_column: str,
    per_reference: int,
    seed: int,
    require_arrow: bool,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any], Sequence[str]]:
    rng = random.Random(seed)

    with input_csv.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"No columns found in CSV: {input_csv}")
        if reference_column not in fieldnames:
            raise ValueError(f"Missing reference column '{reference_column}' in: {input_csv}")
        if reaction_column not in fieldnames:
            raise ValueError(f"Missing reaction column '{reaction_column}' in: {input_csv}")

        grouped: DefaultDict[str, List[Tuple[int, Dict[str, str]]]] = defaultdict(list)
        seen_smiles_per_ref: DefaultDict[str, set[str]] = defaultdict(set)
        total_rows = 0
        valid_rows = 0

        for row_number, row in enumerate(reader, start=2):
            total_rows += 1
            if not isinstance(row, dict):
                continue
            reference = _clean_text(row.get(reference_column))
            if not reference:
                continue

            reaction_smiles = _clean_text(row.get(reaction_column))
            if not _is_valid_reaction_smiles(reaction_smiles, require_arrow=require_arrow):
                continue
            if reaction_smiles in seen_smiles_per_ref[reference]:
                continue

            seen_smiles_per_ref[reference].add(reaction_smiles)
            grouped[reference].append((row_number, row))
            valid_rows += 1

    output_rows: List[Dict[str, Any]] = []
    reference_stats: List[Dict[str, Any]] = []
    references_with_lt_target = 0

    for reference in sorted(grouped.keys()):
        rows = grouped[reference]
        target = max(1, int(per_reference))
        if len(rows) <= target:
            picked = list(rows)
        else:
            picked = rng.sample(rows, k=target)
            picked.sort(key=lambda item: item[0])

        for rank, (row_number, row) in enumerate(picked, start=1):
            payload: Dict[str, Any] = {
                "source_file": input_csv.name,
                "original_row_number": row_number,
                "sample_rank_in_reference": rank,
                "rows_in_reference": len(rows),
            }
            for column in fieldnames:
                payload[column] = row.get(column, "")
            output_rows.append(payload)

        if len(rows) < target:
            references_with_lt_target += 1
        reference_stats.append(
            {
                "reference": reference,
                "rows_in_reference": len(rows),
                "sampled_rows": len(picked),
            }
        )

    output_rows.sort(key=lambda row: (_clean_text(row.get(reference_column)), int(row["sample_rank_in_reference"])))
    summary: Dict[str, Any] = {
        "input_csv": str(input_csv),
        "reference_column": reference_column,
        "reaction_column": reaction_column,
        "seed": seed,
        "per_reference": int(per_reference),
        "total_rows": total_rows,
        "valid_rows": valid_rows,
        "unique_references": len(grouped),
        "references_with_lt_target": references_with_lt_target,
        "sampled_rows_total": len(output_rows),
        "references": reference_stats,
    }
    return output_rows, summary, fieldnames


def _write_csv(path: Path, rows: Iterable[Dict[str, Any]], *, fieldnames: Sequence[str]) -> None:
    output_columns = list(META_COLUMNS) + [c for c in fieldnames if c not in META_COLUMNS]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in output_columns})


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sample deterministic reactions per reference")
    parser.add_argument(
        "--input",
        default="data-processor/reaction_dataset/C_N_Coupling.csv",
        help="Input reaction CSV path",
    )
    parser.add_argument(
        "--output",
        default="results/C_N_Coupling.reference2_test.csv",
        help="Output sampled CSV path",
    )
    parser.add_argument(
        "--summary-output",
        default="results/C_N_Coupling.reference2_test.summary.json",
        help="Output summary JSON path",
    )
    parser.add_argument(
        "--reference-column",
        default="reference",
        help="Reference-group column name",
    )
    parser.add_argument(
        "--reaction-column",
        default="reaction_smiles",
        help="Reaction SMILES column name",
    )
    parser.add_argument(
        "--per-reference",
        type=int,
        default=2,
        help="Samples per reference group",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260211,
        help="Random seed for deterministic sampling",
    )
    parser.add_argument(
        "--allow-non-reaction-smiles",
        action="store_true",
        help="Allow rows without '>>' in reaction_smiles",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_csv = Path(args.input).expanduser().resolve()
    output_csv = Path(args.output).expanduser().resolve()
    summary_json = Path(args.summary_output).expanduser().resolve()

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")
    if int(args.per_reference) <= 0:
        raise ValueError("--per-reference must be > 0")

    rows, summary, fieldnames = build_reference_samples(
        input_csv=input_csv,
        reference_column=str(args.reference_column).strip(),
        reaction_column=str(args.reaction_column).strip(),
        per_reference=int(args.per_reference),
        seed=int(args.seed),
        require_arrow=not bool(args.allow_non_reaction_smiles),
    )
    _write_csv(output_csv, rows, fieldnames=fieldnames)
    _write_json(summary_json, summary)

    print(f"Sampled CSV written: {output_csv}")
    print(f"Summary JSON written: {summary_json}")
    print(f"Unique references: {summary['unique_references']}")
    print(f"Sampled rows: {summary['sampled_rows_total']}")
    print(f"References with < {summary['per_reference']} rows: {summary['references_with_lt_target']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
