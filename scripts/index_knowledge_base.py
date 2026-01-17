#!/usr/bin/env python3
"""
Index the curated knowledge base folder for RAG.

This script reads files under knowledge_base/ and emits a lightweight JSON
summary for downstream retrieval pipelines.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional


_DEFAULT_EXTENSIONS = (".md", ".txt", ".json", ".jsonl", ".yml", ".yaml")


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="utf-8", errors="replace")


def _first_heading(lines: Iterable[str]) -> Optional[str]:
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            return stripped.lstrip("#").strip()
    return None


def _extract_tags(lines: Iterable[str]) -> List[str]:
    tags: List[str] = []
    for line in lines:
        if line.lower().startswith("tags:"):
            tail = line.split(":", 1)[1]
            tags = [t.strip() for t in tail.split(",") if t.strip()]
            break
    return tags


def _index_file(path: Path) -> Dict[str, object]:
    text = _read_text(path)
    lines = text.splitlines()
    heading = _first_heading(lines)
    tags = _extract_tags(lines)
    has_table = any("|" in line for line in lines)
    return {
        "path": str(path.as_posix()),
        "size_bytes": path.stat().st_size,
        "line_count": len(lines),
        "heading": heading or "",
        "tags": tags,
        "has_table": has_table,
    }


def _iter_files(root: Path, extensions: Iterable[str]) -> Iterable[Path]:
    for path in root.rglob("*"):
        if path.is_file() and path.suffix.lower() in extensions:
            yield path


def build_index(root: Path, extensions: Iterable[str]) -> List[Dict[str, object]]:
    return [_index_file(path) for path in sorted(_iter_files(root, extensions))]


def main() -> int:
    parser = argparse.ArgumentParser(description="Index knowledge_base folder.")
    parser.add_argument(
        "--root",
        default="knowledge_base",
        help="Root folder to index (default: knowledge_base).",
    )
    parser.add_argument(
        "--out",
        default="",
        help="Optional output JSON file. When omitted, prints to stdout.",
    )
    parser.add_argument(
        "--extensions",
        default=",".join(_DEFAULT_EXTENSIONS),
        help="Comma-separated file extensions to include.",
    )
    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"Missing knowledge base folder: {root}")

    extensions = tuple(
        ext.strip().lower() for ext in args.extensions.split(",") if ext.strip()
    )
    index = build_index(root, extensions)

    payload = json.dumps(index, indent=2)
    if args.out:
        out_path = Path(args.out)
        out_path.write_text(payload, encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
