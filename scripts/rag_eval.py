#!/usr/bin/env python3
"""
Lightweight RAG evaluation against data/knowledge_base chunks.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Iterable, List


_STOPWORDS = {
    "the", "a", "an", "and", "or", "to", "of", "for", "in", "on", "with",
    "by", "is", "are", "was", "were", "from", "as", "at", "be", "this",
    "that", "it", "its", "into", "about", "top", "list", "notes",
}


def _tokenize(text: str) -> List[str]:
    tokens = re.findall(r"[A-Za-z0-9\-]+", text.lower())
    return [t for t in tokens if t not in _STOPWORDS]


def _score_chunk(query_tokens: List[str], chunk: Dict[str, object]) -> float:
    chunk_tokens = _tokenize(str(chunk.get("text", "")))
    if not chunk_tokens:
        return 0.0
    overlap = len(set(query_tokens) & set(chunk_tokens))
    tag_bonus = 0.0
    tags = chunk.get("tags") or []
    for tag in tags:
        if tag.lower() in query_tokens:
            tag_bonus += 1.5
    return overlap + tag_bonus


def _load_queries(path: Path) -> Iterable[Dict[str, object]]:
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        yield json.loads(line)


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate simple RAG retrieval.")
    parser.add_argument(
        "--root",
        default="data/knowledge_base",
        help="Knowledge base root folder.",
    )
    parser.add_argument(
        "--queries",
        default="data/knowledge_base/eval_queries.jsonl",
        help="Path to JSONL queries file.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=3,
        help="Top-K chunks to consider.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=300,
        help="Approximate words per chunk.",
    )
    parser.add_argument(
        "--chunk-overlap",
        type=int,
        default=40,
        help="Approximate overlap between chunks.",
    )
    parser.add_argument(
        "--exclude",
        default="TEMPLATE.md,eval_queries.jsonl,README.md",
        help="Comma-separated file names to exclude.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in __import__("sys").path:
        __import__("sys").path.insert(0, str(repo_root))
    from scripts.index_knowledge_base import build_index

    index = build_index(
        Path(args.root),
        (".md", ".txt", ".json", ".jsonl", ".yml", ".yaml"),
        chunk_size=max(1, args.chunk_size),
        chunk_overlap=max(0, args.chunk_overlap),
        excludes=[e.strip() for e in args.exclude.split(",") if e.strip()],
    )
    chunks = index["chunks"]

    queries_path = Path(args.queries)
    if not queries_path.exists():
        raise SystemExit(f"Missing queries file: {queries_path}")

    total = 0
    hits = 0
    for entry in _load_queries(queries_path):
        query = str(entry.get("query") or "").strip()
        if not query:
            continue
        total += 1
        expected_doc = str(entry.get("expected_doc") or "")
        expected_tags = [str(t) for t in (entry.get("expected_tags") or [])]
        q_tokens = _tokenize(query)
        scored = [
            (chunk, _score_chunk(q_tokens, chunk)) for chunk in chunks
        ]
        scored.sort(key=lambda item: item[1], reverse=True)
        top_chunks = scored[: max(1, args.top_k)]
        top_paths = [c[0].get("path") for c in top_chunks]
        tag_hits = []
        for chunk, _ in top_chunks:
            tags = [str(t).lower() for t in (chunk.get("tags") or [])]
            for expected in expected_tags:
                if expected.lower() in tags:
                    tag_hits.append(expected)
        hit = (expected_doc in top_paths) or bool(tag_hits)
        if hit:
            hits += 1
        print(f"Query: {query}")
        print(f"- Expected doc: {expected_doc or 'n/a'}")
        print(f"- Top paths: {top_paths}")
        print(f"- Tag hits: {sorted(set(tag_hits))}")
        print(f"- Hit: {hit}")
        print("")

    if total == 0:
        print("No queries evaluated.")
        return 1
    print(f"Hit rate: {hits}/{total} ({(hits / total) * 100:.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
