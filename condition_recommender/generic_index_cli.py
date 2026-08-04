"""CLI for building a deterministic persisted generic reaction index."""

from __future__ import annotations

import argparse
import json

from .generic_indexing import load_generic_index
from .sqlite_indexing import save_sqlite_generic_index


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a versioned generic reaction index from records.jsonl"
    )
    parser.add_argument("records_path", help="Canonical generic records.jsonl")
    parser.add_argument(
        "output_path",
        help="Destination SQLite runtime index (for example generic_index.sqlite)",
    )
    parser.add_argument(
        "--include-review-core",
        action="store_true",
        help=(
            "Build an expert-use trusted-and-review-core index; use a distinct "
            "output such as generic_review_index.sqlite"
        ),
    )
    args = parser.parse_args()
    if not str(args.output_path).casefold().endswith(
        (".sqlite", ".sqlite3", ".db")
    ):
        parser.error(
            "output_path must be a SQLite index; persisted JSON runtime indexes "
            "have been retired"
        )
    index = load_generic_index(
        args.records_path,
        include_review=args.include_review_core,
    )
    report = save_sqlite_generic_index(index, args.output_path)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
