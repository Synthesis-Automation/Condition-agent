"""Read-only migration audit and accepted/quarantine exports."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

from .loader import (
    ADDITIONS_PATH,
    SUBSTANCES_PATH,
    iter_addition_rows,
    iter_substance_rows,
)
from .validation import validate_registry


def _fieldnames(path: Path) -> Tuple[str, ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return tuple(csv.DictReader(handle).fieldnames or ())


def _definition_rows() -> Iterable[Tuple[str, int, Dict[str, str]]]:
    for path, rows in (
        (SUBSTANCES_PATH, iter_substance_rows()),
        (ADDITIONS_PATH, iter_addition_rows()),
    ):
        for row_number, row in enumerate(rows, start=2):
            yield path.name, row_number, row


def run_migration_audit(output_dir: str | Path) -> Dict[str, Any]:
    """Audit both legacy and curated substance definitions without mutating them."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report = validate_registry()
    issues = {
        (str(item.get("definition") or SUBSTANCES_PATH.name), int(item["row_number"])):
        item
        for item in report["issues"]
    }
    union_fields = tuple(
        dict.fromkeys(
            (
                "definition",
                *_fieldnames(SUBSTANCES_PATH),
                *_fieldnames(ADDITIONS_PATH),
            )
        )
    )
    accepted_path = output_dir / "accepted.csv"
    quarantine_path = output_dir / "quarantine.csv"
    with accepted_path.open(
        "w", encoding="utf-8-sig", newline=""
    ) as accepted_handle, quarantine_path.open(
        "w", encoding="utf-8-sig", newline=""
    ) as quarantine_handle:
        accepted_writer = csv.DictWriter(
            accepted_handle,
            fieldnames=union_fields,
        )
        accepted_writer.writeheader()
        quarantine_writer = csv.DictWriter(
            quarantine_handle,
            fieldnames=(*union_fields, "migration_issues"),
        )
        quarantine_writer.writeheader()
        for definition, row_number, row in _definition_rows():
            issue = issues.get((definition, row_number))
            payload = {"definition": definition, **row}
            if issue:
                quarantine_writer.writerow(
                    {
                        **payload,
                        "migration_issues": "|".join(issue["issues"]),
                    }
                )
            else:
                accepted_writer.writerow(payload)
    report_without_rows = {
        key: value
        for key, value in report.items()
        if key not in {"issues", "identifier_issues"}
    }
    (output_dir / "migration_report.json").write_text(
        json.dumps(report_without_rows, indent=2) + "\n",
        encoding="utf-8",
    )
    with (output_dir / "issues.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("definition", "row_number", "name", "cas_raw", "issues"),
        )
        writer.writeheader()
        for item in report["issues"]:
            writer.writerow(
                {
                    **item,
                    "definition": item.get("definition") or SUBSTANCES_PATH.name,
                    "issues": "|".join(item["issues"]),
                }
            )
    return report_without_rows


__all__ = ["run_migration_audit"]
