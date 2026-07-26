"""CLI for exporting a concise reaction-family review CSV."""

from __future__ import annotations

import argparse
import json

from .conversion.concise_review import export_concise_reaction_review_csv


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export five high-value reaction-family review columns"
    )
    parser.add_argument("records_path", help="Canonical records JSONL or JSONL.GZ")
    parser.add_argument("output_path", help="Destination review CSV")
    args = parser.parse_args()
    report = export_concise_reaction_review_csv(
        args.records_path,
        args.output_path,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
