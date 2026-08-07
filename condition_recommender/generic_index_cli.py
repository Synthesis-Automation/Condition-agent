"""CLI for building a deterministic persisted generic reaction index."""

from __future__ import annotations

import argparse
import gzip
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
    if path.name == "shard_manifest.json":
        from .conversion.sharded import iter_gzip_jsonl, validate_sharded_conversion

        integrity = validate_sharded_conversion(path, verify_rows=False)
        if not integrity["valid"]:
            raise ValueError("Sharded conversion integrity check failed")
        manifest = json.loads(path.read_text(encoding="utf-8"))
        for entry in manifest.get("shards") or ():
            yield from iter_gzip_jsonl(path.parent / entry["output_path"])
        return
    opener = gzip.open if path.suffix.casefold() == ".gz" else Path.open
    arguments = (
        {"mode": "rt", "encoding": "utf-8"}
        if path.suffix.casefold() == ".gz"
        else {"mode": "r", "encoding": "utf-8"}
    )
    with opener(path, **arguments) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"JSONL line {line_number} is not an object")
            yield value


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
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
