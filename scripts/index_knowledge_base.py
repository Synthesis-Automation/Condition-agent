#!/usr/bin/env python3
"""
Index the curated knowledge base folder for RAG.

This script reads files under data/knowledge_base/ and emits a lightweight JSON
summary for downstream retrieval pipelines.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


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


def _parse_front_matter(lines: List[str]) -> Tuple[Dict[str, str], int]:
    if not lines or lines[0].strip() != "---":
        return {}, 0
    end_idx = None
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            end_idx = i
            break
    if end_idx is None:
        return {}, 0
    meta: Dict[str, str] = {}
    for line in lines[1:end_idx]:
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        meta[key.strip().lower()] = value.strip()
    return meta, end_idx + 1


def _extract_tags(lines: Iterable[str], meta: Dict[str, str]) -> List[str]:
    tags: List[str] = []
    meta_tags = meta.get("tags", "")
    if meta_tags:
        tags.extend([t.strip() for t in meta_tags.split(",") if t.strip()])
    for line in lines:
        if line.lower().startswith("tags:"):
            tail = line.split(":", 1)[1]
            tags.extend([t.strip() for t in tail.split(",") if t.strip()])
            break
    seen = set()
    out = []
    for tag in tags:
        if tag not in seen:
            out.append(tag)
            seen.add(tag)
    return out


def _split_sections(lines: List[str]) -> List[Tuple[str, str]]:
    sections: List[Tuple[str, str]] = []
    current_heading = ""
    current_lines: List[str] = []
    for line in lines:
        if line.strip().startswith("#"):
            if current_lines:
                sections.append((current_heading, "\n".join(current_lines).strip()))
                current_lines = []
            current_heading = line.lstrip("#").strip()
        current_lines.append(line)
    if current_lines:
        sections.append((current_heading, "\n".join(current_lines).strip()))
    return sections


def _chunk_words(words: List[str], max_words: int, overlap: int) -> List[List[str]]:
    if max_words <= 0:
        return [words]
    chunks: List[List[str]] = []
    start = 0
    while start < len(words):
        end = min(start + max_words, len(words))
        chunks.append(words[start:end])
        if end == len(words):
            break
        start = max(0, end - overlap)
    return chunks


def _index_file(path: Path) -> Dict[str, object]:
    text = _read_text(path)
    lines = text.splitlines()
    meta, body_start = _parse_front_matter(lines)
    body_lines = lines[body_start:] if body_start else lines
    heading = _first_heading(body_lines)
    tags = _extract_tags(lines, meta)
    has_table = any("|" in line for line in lines)
    return {
        "path": str(path.as_posix()),
        "size_bytes": path.stat().st_size,
        "line_count": len(lines),
        "heading": heading or "",
        "tags": tags,
        "has_table": has_table,
        "meta": meta,
    }


def _iter_files(
    root: Path,
    extensions: Iterable[str],
    excludes: Iterable[str],
) -> Iterable[Path]:
    exclude_set = {name.strip() for name in excludes if name.strip()}
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.name in exclude_set:
            continue
        if path.suffix.lower() in extensions:
            yield path


def _build_chunks(
    path: Path,
    *,
    chunk_size: int,
    chunk_overlap: int,
) -> List[Dict[str, object]]:
    text = _read_text(path)
    lines = text.splitlines()
    meta, body_start = _parse_front_matter(lines)
    body_lines = lines[body_start:] if body_start else lines
    tags = _extract_tags(lines, meta)
    heading = _first_heading(body_lines) or ""
    sections = _split_sections(body_lines)
    chunks: List[Dict[str, object]] = []
    for section_heading, section_text in sections:
        words = section_text.split()
        for idx, chunk_words in enumerate(
            _chunk_words(words, max_words=chunk_size, overlap=chunk_overlap),
            start=1,
        ):
            chunk_text = " ".join(chunk_words).strip()
            if not chunk_text:
                continue
            chunk_id = f"{path.as_posix()}::chunk_{len(chunks) + 1}"
            chunks.append(
                {
                    "id": chunk_id,
                    "path": str(path.as_posix()),
                    "heading": heading,
                    "section_heading": section_heading,
                    "tags": tags,
                    "meta": meta,
                    "text": chunk_text,
                    "chunk_index": len(chunks) + 1,
                }
            )
    return chunks


def build_index(
    root: Path,
    extensions: Iterable[str],
    *,
    chunk_size: int,
    chunk_overlap: int,
    excludes: Iterable[str],
) -> Dict[str, object]:
    docs = [
        _index_file(path)
        for path in sorted(_iter_files(root, extensions, excludes))
    ]
    chunks: List[Dict[str, object]] = []
    for doc in docs:
        chunks.extend(
            _build_chunks(
                Path(doc["path"]),
                chunk_size=chunk_size,
                chunk_overlap=chunk_overlap,
            )
        )
    return {"docs": docs, "chunks": chunks}


def main() -> int:
    parser = argparse.ArgumentParser(description="Index data/knowledge_base folder.")
    parser.add_argument(
        "--root",
        default="data/knowledge_base",
        help="Root folder to index (default: data/knowledge_base).",
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
    parser.add_argument(
        "--exclude",
        default="TEMPLATE.md,eval_queries.jsonl,README.md",
        help="Comma-separated file names to exclude.",
    )
    parser.add_argument(
        "--mode",
        default="both",
        choices=("summary", "chunks", "both"),
        help="Output summary, chunks, or both (default: both).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=300,
        help="Approximate words per chunk (default: 300).",
    )
    parser.add_argument(
        "--chunk-overlap",
        type=int,
        default=40,
        help="Approximate word overlap between chunks (default: 40).",
    )
    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"Missing knowledge base folder: {root}")

    extensions = tuple(
        ext.strip().lower() for ext in args.extensions.split(",") if ext.strip()
    )
    index = build_index(
        root,
        extensions,
        chunk_size=max(1, args.chunk_size),
        chunk_overlap=max(0, args.chunk_overlap),
        excludes=[e.strip() for e in args.exclude.split(",") if e.strip()],
    )
    if args.mode == "summary":
        payload_obj = index["docs"]
    elif args.mode == "chunks":
        payload_obj = index["chunks"]
    else:
        payload_obj = index

    payload = json.dumps(payload_obj, indent=2)
    if args.out:
        out_path = Path(args.out)
        out_path.write_text(payload, encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
