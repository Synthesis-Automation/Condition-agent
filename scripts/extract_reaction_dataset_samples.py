#!/usr/bin/env python3
"""
Extract deterministic reaction samples from each CSV in reaction_dataset.

Example:
  python scripts/extract_reaction_dataset_samples.py \
    --input-dir data-processor/reaction_dataset \
    --output examples/test_reactions_from_reaction_dataset.csv \
    --summary-output results/reaction_dataset_sample_summary.json \
    --per-file 3
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


OUTPUT_COLUMNS = [
    "source_file",
    "source_reaction_family",
    "sample_rank_in_file",
    "original_row_number",
    "reaction_id",
    "reaction_type",
    "reaction_smiles",
    "stages",
    "steps",
    "reference",
    "notes",
]


def _to_int(value: Any, default: int = 0) -> int:
    try:
        text = str(value or "").strip()
        if not text:
            return default
        return int(float(text))
    except Exception:
        return default


def _is_multistep(row: Dict[str, str]) -> bool:
    stages = _to_int(row.get("stages"), 0)
    steps = _to_int(row.get("steps"), 0)
    if stages > 1 or steps > 1:
        return True
    notes = str(row.get("notes") or "").strip().lower()
    return "step" in notes or "multi" in notes


def _is_valid_reaction_smiles(text: str, require_arrow: bool) -> bool:
    value = str(text or "").strip()
    if not value:
        return False
    if require_arrow and ">>" not in value:
        return False
    return True


def _pick_rows_for_file(
    rows: List[Tuple[int, Dict[str, str]]],
    *,
    per_file: int,
    rng: random.Random,
) -> List[Tuple[int, Dict[str, str]]]:
    if len(rows) <= per_file:
        return list(rows)

    multi = [row for row in rows if _is_multistep(row[1])]
    non_multi = [row for row in rows if not _is_multistep(row[1])]

    selected: List[Tuple[int, Dict[str, str]]] = []
    if multi:
        selected.append(rng.choice(multi))

    needed = per_file - len(selected)
    remaining = [row for row in rows if row not in selected]
    if needed > 0:
        selected.extend(rng.sample(remaining, k=needed))

    selected.sort(key=lambda item: item[0])
    return selected


def _iter_dataset_csv_files(input_dir: Path) -> List[Path]:
    return sorted(path for path in input_dir.glob("*.csv") if path.is_file())


def build_samples(
    *,
    input_dir: Path,
    per_file: int,
    seed: int,
    require_arrow: bool,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    rng = random.Random(seed)
    output_rows: List[Dict[str, Any]] = []

    files = _iter_dataset_csv_files(input_dir)
    summary_files: List[Dict[str, Any]] = []
    total_valid_rows = 0

    for csv_path in files:
        with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "reaction_smiles" not in reader.fieldnames:
                summary_files.append(
                    {
                        "file": csv_path.name,
                        "status": "skipped_missing_reaction_smiles",
                        "total_rows": 0,
                        "valid_rows": 0,
                        "sampled_rows": 0,
                    }
                )
                continue

            all_valid_rows: List[Tuple[int, Dict[str, str]]] = []
            seen_smiles: set[str] = set()
            for row_index, row in enumerate(reader, start=2):
                rxn = str((row or {}).get("reaction_smiles") or "").strip()
                if not _is_valid_reaction_smiles(rxn, require_arrow=require_arrow):
                    continue
                if rxn in seen_smiles:
                    continue
                seen_smiles.add(rxn)
                all_valid_rows.append((row_index, row))

            picked = _pick_rows_for_file(all_valid_rows, per_file=per_file, rng=rng)
            total_valid_rows += len(all_valid_rows)

            family = csv_path.stem
            for rank, (row_index, row) in enumerate(picked, start=1):
                output_rows.append(
                    {
                        "source_file": csv_path.name,
                        "source_reaction_family": family,
                        "sample_rank_in_file": rank,
                        "original_row_number": row_index,
                        "reaction_id": str(row.get("reaction_id") or "").strip(),
                        "reaction_type": str(row.get("reaction_type") or "").strip(),
                        "reaction_smiles": str(row.get("reaction_smiles") or "").strip(),
                        "stages": str(row.get("stages") or "").strip(),
                        "steps": str(row.get("steps") or "").strip(),
                        "reference": str(row.get("reference") or "").strip(),
                        "notes": str(row.get("notes") or "").strip(),
                    }
                )

            summary_files.append(
                {
                    "file": csv_path.name,
                    "status": "ok",
                    "total_rows": len(all_valid_rows),
                    "valid_rows": len(all_valid_rows),
                    "sampled_rows": len(picked),
                    "multistep_in_sample": sum(1 for _, row in picked if _is_multistep(row)),
                }
            )

    output_rows.sort(key=lambda row: (row["source_file"], int(row["sample_rank_in_file"])))
    summary = {
        "input_dir": str(input_dir),
        "seed": seed,
        "per_file": per_file,
        "files_scanned": len(files),
        "total_samples": len(output_rows),
        "total_valid_rows": total_valid_rows,
        "files": summary_files,
    }
    return output_rows, summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract sampled reactions per dataset file")
    parser.add_argument(
        "--input-dir",
        default="data-processor/reaction_dataset",
        help="Directory with reaction CSV files",
    )
    parser.add_argument(
        "--output",
        default="examples/test_reactions_from_reaction_dataset.csv",
        help="Output sampled CSV path",
    )
    parser.add_argument(
        "--summary-output",
        default="results/reaction_dataset_sample_summary.json",
        help="Output summary JSON path",
    )
    parser.add_argument(
        "--per-file",
        type=int,
        default=3,
        help="Number of samples per input CSV file",
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


def _write_csv(path: Path, rows: Iterable[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in OUTPUT_COLUMNS})


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main() -> int:
    args = _parse_args()
    input_dir = Path(args.input_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    summary_output_path = Path(args.summary_output).expanduser().resolve()

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if args.per_file <= 0:
        raise ValueError("--per-file must be > 0")

    rows, summary = build_samples(
        input_dir=input_dir,
        per_file=int(args.per_file),
        seed=int(args.seed),
        require_arrow=not bool(args.allow_non_reaction_smiles),
    )
    _write_csv(output_path, rows)
    _write_json(summary_output_path, summary)

    print(f"Sampled CSV written: {output_path}")
    print(f"Summary JSON written: {summary_output_path}")
    print(f"Files scanned: {summary['files_scanned']}")
    print(f"Total samples: {summary['total_samples']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
