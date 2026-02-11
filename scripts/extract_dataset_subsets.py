#!/usr/bin/env python3
"""
Extract deterministic per-dataset subsets into a separate folder.

Example:
  python scripts/extract_dataset_subsets.py \
    --input-dir data-processor/reaction_dataset \
    --output-dir results/reaction_dataset_300 \
    --per-file 300
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import zlib
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


def _is_valid_reaction_smiles(value: Any, *, require_arrow: bool) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    if require_arrow and ">>" not in text:
        return False
    return True


def _iter_csv_files(input_dir: Path) -> List[Path]:
    return sorted(path for path in input_dir.glob("*.csv") if path.is_file())


def _pick_rows(
    rows: List[Tuple[int, Dict[str, str]]],
    *,
    per_file: int,
    seed: int,
    file_name: str,
) -> List[Tuple[int, Dict[str, str]]]:
    if len(rows) <= per_file:
        return list(rows)
    file_seed = int(seed) + int(zlib.crc32(file_name.encode("utf-8")))
    rng = random.Random(file_seed)
    picked = rng.sample(rows, k=per_file)
    picked.sort(key=lambda item: item[0])
    return picked


def _write_csv(
    path: Path,
    *,
    fieldnames: Iterable[str],
    rows: Iterable[Dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_dataset_subsets(
    *,
    input_dir: Path,
    output_dir: Path,
    per_file: int,
    seed: int,
    reaction_column: str,
    require_arrow: bool,
) -> Dict[str, Any]:
    files = _iter_csv_files(input_dir)
    summary_files: List[Dict[str, Any]] = []

    total_input_rows = 0
    total_valid_rows = 0
    total_sampled_rows = 0
    files_with_lt_target = 0

    for csv_path in files:
        with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or [])
            if not fieldnames:
                summary_files.append(
                    {
                        "file": csv_path.name,
                        "status": "skipped_empty_header",
                        "input_rows": 0,
                        "valid_rows": 0,
                        "sampled_rows": 0,
                    }
                )
                continue
            if reaction_column not in fieldnames:
                summary_files.append(
                    {
                        "file": csv_path.name,
                        "status": f"skipped_missing_{reaction_column}",
                        "input_rows": 0,
                        "valid_rows": 0,
                        "sampled_rows": 0,
                    }
                )
                continue

            input_rows = 0
            valid_rows: List[Tuple[int, Dict[str, str]]] = []
            seen_smiles: set[str] = set()
            for row_number, row in enumerate(reader, start=2):
                input_rows += 1
                if not isinstance(row, dict):
                    continue
                rxn = str(row.get(reaction_column) or "").strip()
                if not _is_valid_reaction_smiles(rxn, require_arrow=require_arrow):
                    continue
                if rxn in seen_smiles:
                    continue
                seen_smiles.add(rxn)
                valid_rows.append((row_number, row))

            picked = _pick_rows(
                valid_rows,
                per_file=max(1, int(per_file)),
                seed=int(seed),
                file_name=csv_path.name,
            )
            sampled_rows = [row for _, row in picked]

            out_path = output_dir / csv_path.name
            _write_csv(out_path, fieldnames=fieldnames, rows=sampled_rows)

            if len(valid_rows) < max(1, int(per_file)):
                files_with_lt_target += 1

            total_input_rows += input_rows
            total_valid_rows += len(valid_rows)
            total_sampled_rows += len(sampled_rows)
            summary_files.append(
                {
                    "file": csv_path.name,
                    "status": "ok",
                    "input_rows": input_rows,
                    "valid_rows": len(valid_rows),
                    "sampled_rows": len(sampled_rows),
                }
            )

    summary_files.sort(key=lambda item: str(item.get("file") or ""))
    summary = {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "per_file": int(per_file),
        "seed": int(seed),
        "reaction_column": reaction_column,
        "files_scanned": len(files),
        "files_with_lt_target": files_with_lt_target,
        "total_input_rows": total_input_rows,
        "total_valid_rows": total_valid_rows,
        "total_sampled_rows": total_sampled_rows,
        "files": summary_files,
    }
    return summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract per-dataset subsets into a folder")
    parser.add_argument(
        "--input-dir",
        default="data-processor/reaction_dataset",
        help="Directory containing reaction CSV files",
    )
    parser.add_argument(
        "--output-dir",
        default="results/reaction_dataset_300",
        help="Output directory for per-file sampled CSVs",
    )
    parser.add_argument(
        "--summary-output",
        default="results/reaction_dataset_300_summary.json",
        help="Summary JSON output path",
    )
    parser.add_argument(
        "--per-file",
        type=int,
        default=300,
        help="Rows sampled per dataset file",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260211,
        help="Random seed for deterministic sampling",
    )
    parser.add_argument(
        "--reaction-column",
        default="reaction_smiles",
        help="Reaction smiles column name in input CSVs",
    )
    parser.add_argument(
        "--allow-non-reaction-smiles",
        action="store_true",
        help="Allow rows without '>>' in reaction smiles",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    summary_output = Path(args.summary_output).expanduser().resolve()

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if int(args.per_file) <= 0:
        raise ValueError("--per-file must be > 0")

    summary = build_dataset_subsets(
        input_dir=input_dir,
        output_dir=output_dir,
        per_file=int(args.per_file),
        seed=int(args.seed),
        reaction_column=str(args.reaction_column).strip(),
        require_arrow=not bool(args.allow_non_reaction_smiles),
    )
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(f"Output folder: {output_dir}")
    print(f"Summary JSON: {summary_output}")
    print(f"Files scanned: {summary['files_scanned']}")
    print(f"Total sampled rows: {summary['total_sampled_rows']}")
    print(f"Files with < {summary['per_file']} valid rows: {summary['files_with_lt_target']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
