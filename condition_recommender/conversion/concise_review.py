"""Concise chemistry review export for canonical generic conversion records."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, Mapping, Optional

from .generic import GenericConversionCache, convert_record
from .input_schema import discover_csv_datasets, iter_csv_records

CONCISE_REACTION_REVIEW_FIELDS = (
    "canonical_reaction_smiles",
    "reaction_display_label_detailed",
    "original_reaction_type",
    "detected_reaction_family",
    "detection_status",
)


@dataclass(frozen=True)
class ConciseReviewProgress:
    """Progress update for a recursive dataset-folder review export."""

    phase: str
    file_index: int
    file_count: int
    row_count: int
    current_file: str
    message: str


class ConciseReviewConversionCancelled(RuntimeError):
    """Raised when a caller cancels a folder conversion between records."""


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


def convert_dataset_folder_to_concise_review_csv(
    dataset_path: str | Path,
    output_path: str | Path,
    *,
    progress_callback: Optional[Callable[[ConciseReviewProgress], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
    progress_interval: int = 25,
) -> Dict[str, Any]:
    """Recursively convert source CSVs directly into the concise review view."""
    if progress_interval < 1:
        raise ValueError("progress_interval must be positive")
    source = Path(dataset_path)
    destination = Path(output_path)
    paths = tuple(
        path
        for path in discover_csv_datasets(source)
        if path.resolve() != destination.resolve()
    )
    if not paths:
        raise ValueError(f"No CSV datasets found at {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    row_count = 0
    cache = GenericConversionCache(max_entries=5_000)

    def report(
        phase: str,
        *,
        file_index: int,
        current_file: str,
        message: str,
    ) -> None:
        if progress_callback is not None:
            progress_callback(
                ConciseReviewProgress(
                    phase=phase,
                    file_index=file_index,
                    file_count=len(paths),
                    row_count=row_count,
                    current_file=current_file,
                    message=message,
                )
            )

    report(
        "discovered",
        file_index=0,
        current_file="",
        message=f"Found {len(paths)} CSV dataset file(s).",
    )
    completed = False
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=CONCISE_REACTION_REVIEW_FIELDS,
            )
            writer.writeheader()
            for file_index, path in enumerate(paths, start=1):
                if cancel_check is not None and cancel_check():
                    raise ConciseReviewConversionCancelled(
                        "Review conversion cancelled"
                    )
                relative_path = (
                    path.relative_to(source).as_posix()
                    if source.is_dir()
                    else path.name
                )
                report(
                    "file_started",
                    file_index=file_index,
                    current_file=relative_path,
                    message=(
                        f"Processing file {file_index}/{len(paths)}: "
                        f"{relative_path}"
                    ),
                )
                file_row_count = 0
                for raw_record in iter_csv_records(path):
                    if cancel_check is not None and cancel_check():
                        raise ConciseReviewConversionCancelled(
                            "Review conversion cancelled"
                        )
                    converted = convert_record(raw_record, cache=cache)
                    writer.writerow(
                        concise_reaction_review_row(converted.to_dict())
                    )
                    row_count += 1
                    file_row_count += 1
                    if file_row_count % progress_interval == 0:
                        report(
                            "rows_converted",
                            file_index=file_index,
                            current_file=relative_path,
                            message=(
                                f"{relative_path}: {file_row_count} row(s); "
                                f"{row_count} total."
                            ),
                        )
                report(
                    "file_completed",
                    file_index=file_index,
                    current_file=relative_path,
                    message=(
                        f"Completed {relative_path}: "
                        f"{file_row_count} reaction(s)."
                    ),
                )
        temporary.replace(destination)
        completed = True
    finally:
        if not completed and temporary.is_file():
            temporary.unlink()
    report(
        "completed",
        file_index=len(paths),
        current_file="",
        message=f"Finished {row_count} reaction(s).",
    )
    return {
        "schema_version": "1.0",
        "artifact_type": "concise_reaction_review_csv",
        "dataset_path": str(source),
        "output_path": str(destination),
        "source_file_count": len(paths),
        "row_count": row_count,
        "columns": list(CONCISE_REACTION_REVIEW_FIELDS),
    }


__all__ = [
    "CONCISE_REACTION_REVIEW_FIELDS",
    "ConciseReviewConversionCancelled",
    "ConciseReviewProgress",
    "concise_reaction_review_row",
    "convert_dataset_folder_to_concise_review_csv",
    "export_concise_reaction_review_csv",
]
