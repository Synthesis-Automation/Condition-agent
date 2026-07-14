"""CLI for mixed reaction-dataset quality and chemistry auditing."""

from __future__ import annotations

import argparse
import json

from .conversion.audit import audit_datasets


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit mixed reaction CSVs without modifying source files"
    )
    parser.add_argument("dataset_path", help="CSV file or directory of CSV files")
    parser.add_argument("output_dir", help="Directory for JSON and Markdown reports")
    parser.add_argument(
        "--chemistry-sample-per-file",
        type=int,
        default=100,
        help="Deterministic rows per file for full reaction featurization",
    )
    args = parser.parse_args()
    report = audit_datasets(
        args.dataset_path,
        args.output_dir,
        chemistry_sample_per_file=max(0, args.chemistry_sample_per_file),
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
