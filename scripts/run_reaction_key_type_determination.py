#!/usr/bin/env python3
"""
Run reaction-key and reaction-type determination over a CSV dataset.

Example:
  python scripts/run_reaction_key_type_determination.py \
    --input results/C_N_Coupling.reference2_test.csv \
    --output-csv results/C_N_Coupling.reference2_test.predicted.csv \
    --summary-json results/C_N_Coupling.reference2_test.predicted.summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction


ADDED_COLUMNS = [
    "predicted_reaction_type",
    "predicted_reaction_confidence",
    "predicted_reaction_key",
    "predicted_reaction_key_quality_level",
    "predicted_reaction_key_quality_score",
    "predicted_mapping_warning",
    "prediction_error",
]


def _reaction_type_label(payload: Dict[str, Any]) -> Tuple[str, float]:
    value = payload.get("reaction_type")
    if isinstance(value, dict):
        label = str(value.get("reaction_type") or value.get("name") or "Unknown").strip() or "Unknown"
        try:
            confidence = float(value.get("confidence") or 0.0)
        except Exception:
            confidence = 0.0
        return label, confidence
    text = str(value or "").strip()
    return (text if text else "Unknown"), 0.0


def _extract_key_quality(payload: Dict[str, Any]) -> Tuple[str, float]:
    events = payload.get("reaction_events")
    quality: Dict[str, Any] = {}
    if isinstance(events, dict):
        candidate = events.get("reaction_key_quality")
        if isinstance(candidate, dict):
            quality = candidate
    if not quality:
        detection = payload.get("detection")
        if isinstance(detection, dict):
            candidate = detection.get("reaction_key_quality")
            if isinstance(candidate, dict):
                quality = candidate
    level = str(quality.get("level") or "unknown").strip().lower()
    try:
        score = float(quality.get("score_0_1") or 0.0)
    except Exception:
        score = 0.0
    return level or "unknown", score


def _mapping_warning(payload: Dict[str, Any]) -> bool:
    detection = payload.get("detection")
    if not isinstance(detection, dict):
        return False
    warning = detection.get("mapping_warning")
    return isinstance(warning, dict) and bool(warning)


def run_prediction(
    *,
    input_csv: Path,
    output_csv: Path,
    summary_json: Path,
    reaction_column: str,
    progress_every: int,
) -> Dict[str, Any]:
    with input_csv.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"No columns found in CSV: {input_csv}")
        if reaction_column not in fieldnames:
            raise ValueError(f"Missing reaction column '{reaction_column}' in: {input_csv}")

        output_csv.parent.mkdir(parents=True, exist_ok=True)
        output_columns = list(fieldnames) + [c for c in ADDED_COLUMNS if c not in fieldnames]
        rows_out: List[Dict[str, Any]] = []

        total_rows = 0
        processed = 0
        errors = 0
        unknown = 0
        empty_key = 0
        low_quality = 0
        mapping_warning_count = 0

        options = {
            "detailed": False,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        }

        for row in reader:
            total_rows += 1
            if not isinstance(row, dict):
                continue

            output_row = dict(row)
            for col in ADDED_COLUMNS:
                output_row[col] = ""

            reaction_smiles = str(row.get(reaction_column) or "").strip()
            if not reaction_smiles:
                output_row["prediction_error"] = "missing_reaction_smiles"
                errors += 1
                rows_out.append(output_row)
                continue

            try:
                result = featurize_reaction(reaction_smiles, options=options)
            except Exception as exc:
                output_row["prediction_error"] = f"featurize_error:{type(exc).__name__}"
                errors += 1
                rows_out.append(output_row)
                continue

            processed += 1
            reaction_type, confidence = _reaction_type_label(result)
            reaction_key = str(result.get("reaction_key") or "").strip()
            quality_level, quality_score = _extract_key_quality(result)
            has_mapping_warning = _mapping_warning(result)

            output_row["predicted_reaction_type"] = reaction_type
            output_row["predicted_reaction_confidence"] = round(confidence, 4)
            output_row["predicted_reaction_key"] = reaction_key
            output_row["predicted_reaction_key_quality_level"] = quality_level
            output_row["predicted_reaction_key_quality_score"] = round(quality_score, 4)
            output_row["predicted_mapping_warning"] = "1" if has_mapping_warning else "0"

            if reaction_type.lower() == "unknown":
                unknown += 1
            if not reaction_key:
                empty_key += 1
            if quality_level == "low" or quality_score < 0.45:
                low_quality += 1
            if has_mapping_warning:
                mapping_warning_count += 1

            rows_out.append(output_row)
            if progress_every > 0 and processed % progress_every == 0:
                print(
                    f"[progress] processed={processed} unknown={unknown} "
                    f"empty_key={empty_key} low_quality={low_quality} errors={errors}"
                )

    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_columns)
        writer.writeheader()
        for row in rows_out:
            writer.writerow({key: row.get(key, "") for key in output_columns})

    denom = max(processed, 1)
    summary = {
        "input_csv": str(input_csv),
        "output_csv": str(output_csv),
        "reaction_column": reaction_column,
        "total_rows": total_rows,
        "processed_rows": processed,
        "prediction_errors": errors,
        "unknown_reaction_type": {"count": unknown, "rate": round(unknown / denom, 4)},
        "empty_reaction_key": {"count": empty_key, "rate": round(empty_key / denom, 4)},
        "low_reaction_key_quality": {"count": low_quality, "rate": round(low_quality / denom, 4)},
        "with_mapping_warning": {"count": mapping_warning_count, "rate": round(mapping_warning_count / denom, 4)},
    }

    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run reaction key/type determination on a CSV dataset")
    parser.add_argument(
        "--input",
        default="results/C_N_Coupling.reference2_test.csv",
        help="Input CSV path",
    )
    parser.add_argument(
        "--output-csv",
        default="results/C_N_Coupling.reference2_test.predicted.csv",
        help="Output CSV path with predictions",
    )
    parser.add_argument(
        "--summary-json",
        default="results/C_N_Coupling.reference2_test.predicted.summary.json",
        help="Summary JSON output path",
    )
    parser.add_argument(
        "--reaction-column",
        default="reaction_smiles",
        help="Reaction smiles column in input CSV",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
        help="Print progress every N processed rows (0 to disable)",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_csv = Path(args.input).expanduser().resolve()
    output_csv = Path(args.output_csv).expanduser().resolve()
    summary_json = Path(args.summary_json).expanduser().resolve()

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    summary = run_prediction(
        input_csv=input_csv,
        output_csv=output_csv,
        summary_json=summary_json,
        reaction_column=str(args.reaction_column).strip(),
        progress_every=max(0, int(args.progress_every)),
    )
    print(f"Predicted CSV: {output_csv}")
    print(f"Summary JSON: {summary_json}")
    print(f"Processed rows: {summary['processed_rows']}")
    print(f"Prediction errors: {summary['prediction_errors']}")
    print(f"Unknown reaction types: {summary['unknown_reaction_type']['count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
