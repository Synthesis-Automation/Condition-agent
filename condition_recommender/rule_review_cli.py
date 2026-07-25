"""CLI for exporting a flat chemist-review view of expert condition rules."""

from __future__ import annotations

import argparse

from .rule_review import export_rule_review_csv


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export validated expert rules and explicit condition variants "
            "as a flat review CSV"
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    export_parser = subparsers.add_parser(
        "export",
        help="Write the current validated rule definitions as CSV",
    )
    export_parser.add_argument("output", help="Destination CSV path")
    args = parser.parse_args()
    if args.command == "export":
        rows = export_rule_review_csv(args.output)
        column_count = len(rows[0]) if rows else 0
        print(
            f"Wrote {len(rows)} rule-variant rows with {column_count} "
            f"columns to {args.output}"
        )


if __name__ == "__main__":
    main()
