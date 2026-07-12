"""Read-only migration audit and accepted/quarantine exports."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict

from .loader import SUBSTANCES_PATH, iter_substance_rows
from .validation import validate_registry


def run_migration_audit(output_dir: str | Path) -> Dict[str, Any]:
    output_dir = Path(output_dir); output_dir.mkdir(parents=True, exist_ok=True)
    report = validate_registry(); issues = {int(item["row_number"]): item for item in report["issues"]}
    with SUBSTANCES_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        fieldnames = list(csv.DictReader(handle).fieldnames or [])
    accepted_path, quarantine_path = output_dir / "accepted.csv", output_dir / "quarantine.csv"
    with accepted_path.open("w", encoding="utf-8-sig", newline="") as accepted_handle, quarantine_path.open("w", encoding="utf-8-sig", newline="") as quarantine_handle:
        accepted_writer = csv.DictWriter(accepted_handle, fieldnames=fieldnames); accepted_writer.writeheader()
        quarantine_writer = csv.DictWriter(quarantine_handle, fieldnames=[*fieldnames, "migration_issues"]); quarantine_writer.writeheader()
        for row_number, row in enumerate(iter_substance_rows(), start=2):
            issue = issues.get(row_number)
            if issue:
                quarantine_writer.writerow({**row, "migration_issues": "|".join(issue["issues"])})
            else:
                accepted_writer.writerow(row)
    report_without_rows = {key: value for key, value in report.items() if key != "issues"}
    (output_dir / "migration_report.json").write_text(json.dumps(report_without_rows, indent=2) + "\n", encoding="utf-8")
    with (output_dir / "issues.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("row_number", "name", "cas_raw", "issues")); writer.writeheader()
        for item in report["issues"]: writer.writerow({**item, "issues": "|".join(item["issues"])})
    return report_without_rows


__all__ = ["run_migration_audit"]
