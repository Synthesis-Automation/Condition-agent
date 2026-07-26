"""Concise chemistry review export for canonical generic conversion records."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping

CONCISE_REACTION_REVIEW_FIELDS = (
    "canonical_reaction_smiles",
    "reaction_display_label_detailed",
    "original_reaction_type",
    "detected_reaction_family",
    "detection_status",
)


def _iter_records(path: Path) -> Iterator[Dict[str, Any]]:
    handle = (
        gzip.open(path, mode="rt", encoding="utf-8")
        if path.suffix.casefold() == ".gz"
        else path.open(mode="r", encoding="utf-8")
    )
    with handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid record JSONL at {path}:{line_number}: {exc.msg}"
                ) from exc
            if not isinstance(value, dict):
                raise ValueError(
                    f"Record is not a JSON object at {path}:{line_number}"
                )
            yield value


def concise_reaction_review_row(record: Mapping[str, Any]) -> Dict[str, str]:
    """Select the five fields needed for rapid reaction-family review."""
    display = record.get("reaction_display_label")
    display_value = display if isinstance(display, Mapping) else {}
    return {
        "canonical_reaction_smiles": str(
            record.get("canonical_reaction_smiles")
            or record.get("reaction_smiles")
            or ""
        ),
        "reaction_display_label_detailed": str(
            display_value.get("detailed") or ""
        ),
        "original_reaction_type": str(
            record.get("source_declared_family") or ""
        ),
        "detected_reaction_family": str(record.get("named_family") or ""),
        "detection_status": str(
            display_value.get("status")
            or record.get("reaction_label_status")
            or "unavailable"
        ),
    }


def export_concise_reaction_review_csv(
    records_path: str | Path,
    output_path: str | Path,
) -> Dict[str, Any]:
    """Write an Excel-friendly concise review CSV from canonical JSONL records."""
    source = Path(records_path)
    if not source.is_file():
        raise ValueError(f"Canonical records file does not exist: {source}")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with destination.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CONCISE_REACTION_REVIEW_FIELDS,
        )
        writer.writeheader()
        for record in _iter_records(source):
            writer.writerow(concise_reaction_review_row(record))
            row_count += 1
    return {
        "schema_version": "1.0",
        "artifact_type": "concise_reaction_review_csv",
        "records_path": str(source),
        "output_path": str(destination),
        "row_count": row_count,
        "columns": list(CONCISE_REACTION_REVIEW_FIELDS),
    }


__all__ = [
    "CONCISE_REACTION_REVIEW_FIELDS",
    "concise_reaction_review_row",
    "export_concise_reaction_review_csv",
]
