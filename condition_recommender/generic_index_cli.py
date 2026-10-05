"""CLI for building a deterministic persisted generic reaction index."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable

from .generic_indexing import load_generic_index
from .sqlite_indexing import (
    build_sqlite_generic_index,
    save_sqlite_generic_index,
)


def _iter_records(path: Path) -> Iterable[Dict[str, Any]]:
    """Stream canonical JSONL records without materializing the corpus."""
    from .corpus_io import iter_canonical_records

    yield from iter_canonical_records(path, strict=True)



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
        "--baseline-only",
        action="store_true",
        help="Skip the default shared-core companion for baseline evaluation",
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
    source = Path(args.records_path)
    if source.suffix.casefold() in {".sqlite", ".sqlite3", ".db"}:
        index = load_generic_index(
            source,
            include_review=args.include_review_core,
        )
        report = save_sqlite_generic_index(index, args.output_path)
    else:
        report = build_sqlite_generic_index(
            _iter_records(source),
            args.output_path,
            include_review=args.include_review_core,
        )
    if not args.baseline_only:
        from .shared_core_index import build_shared_core_index

        index = load_generic_index(
            args.output_path, include_review=args.include_review_core
        )
        report["shared_core_report"] = build_shared_core_index(
            index, Path(args.output_path).with_suffix(".shared_core.sqlite")
        )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
