"""Read-only audit and review export for unified substance records."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict

from .loader import SUBSTANCES_PATH, iter_substance_records
from .validation import validate_registry


def run_migration_audit(output_dir: str | Path) -> Dict[str, Any]:
    """Write accepted/quarantined JSONL plus a compact issue CSV."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    report = validate_registry()
    issue_lines = {int(item["line_number"]): item for item in report["issues"]}
    accepted_path = output / "accepted.jsonl"
    quarantine_path = output / "quarantine.jsonl"
    with accepted_path.open("w", encoding="utf-8", newline="\n") as accepted, quarantine_path.open("w", encoding="utf-8", newline="\n") as quarantined:
        for line_number, record in iter_substance_records(SUBSTANCES_PATH):
            target = quarantined if line_number in issue_lines else accepted
            target.write(json.dumps(record, ensure_ascii=True, sort_keys=True) + "\n")
    with (output / "issues.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("line_number", "substance_id", "issues"),
        )
        writer.writeheader()
        for issue in report["issues"]:
            writer.writerow(
                {
                    "line_number": issue.get("line_number", ""),
                    "substance_id": issue.get("substance_id", ""),
                    "issues": "|".join(issue.get("issues") or ()),
                }
            )
    (output / "migration_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = ["run_migration_audit"]
