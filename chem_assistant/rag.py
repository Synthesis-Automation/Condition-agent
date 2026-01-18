"""
Lightweight RAG utilities over the curated data/knowledge_base folder.
"""

from __future__ import annotations

from functools import lru_cache
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


_STOPWORDS = {
    "the",
    "a",
    "an",
    "and",
    "or",
    "to",
    "of",
    "for",
    "in",
    "on",
    "with",
    "by",
    "is",
    "are",
    "was",
    "were",
    "from",
    "as",
    "at",
    "be",
    "this",
    "that",
    "it",
    "its",
    "into",
    "about",
    "top",
    "list",
    "notes",
}


def _tokenize(text: str) -> List[str]:
    tokens = re.findall(r"[A-Za-z0-9\\-]+", text.lower())
    return [t for t in tokens if t not in _STOPWORDS]


def _collect_tags(chunk: Dict[str, object]) -> List[str]:
    tags: List[str] = []
    raw_tags = chunk.get("tags") or []
    if isinstance(raw_tags, list):
        tags.extend([str(t) for t in raw_tags])
    meta = chunk.get("meta") or {}
    if isinstance(meta, dict):
        meta_tags = meta.get("tags")
        if isinstance(meta_tags, str):
            tags.extend([t.strip() for t in meta_tags.split(",") if t.strip()])
    return tags


def _score_chunk(query_tokens: List[str], chunk: Dict[str, object]) -> float:
    chunk_tokens = _tokenize(str(chunk.get("text", "")))
    if not chunk_tokens:
        return 0.0
    overlap = len(set(query_tokens) & set(chunk_tokens))
    tag_bonus = 0.0
    for tag in _collect_tags(chunk):
        if tag.lower() in query_tokens:
            tag_bonus += 1.5
    meta = chunk.get("meta") or {}
    if isinstance(meta, dict):
        doc_id = meta.get("id")
        if isinstance(doc_id, str) and doc_id.lower() in query_tokens:
            tag_bonus += 2.0
    return overlap + tag_bonus


@lru_cache(maxsize=8)
def _load_index(
    root: str,
    chunk_size: int,
    chunk_overlap: int,
    excludes: Tuple[str, ...],
) -> Dict[str, object]:
    from scripts.index_knowledge_base import build_index

    root_path = Path(root)
    extensions = (".md", ".txt", ".json", ".jsonl", ".yml", ".yaml")
    return build_index(
        root_path,
        extensions,
        chunk_size=chunk_size,
        chunk_overlap=chunk_overlap,
        excludes=excludes,
    )


def search_knowledge_base(
    query: str,
    *,
    top_k: int = 3,
    root: Optional[str] = None,
    chunk_size: int = 300,
    chunk_overlap: int = 40,
    include_text: bool = True,
    max_chars: int = 1200,
    excludes: Optional[Iterable[str]] = None,
) -> Dict[str, object]:
    if not query or not query.strip():
        return {"query": query, "results": []}

    repo_root = Path(__file__).resolve().parents[1]
    kb_root = Path(root) if root else (repo_root / "data" / "knowledge_base")
    if not kb_root.exists():
        raise FileNotFoundError(f"Knowledge base not found: {kb_root}")

    exclude_list = tuple(excludes) if excludes is not None else (
        "TEMPLATE.md",
        "eval_queries.jsonl",
        "README.md",
    )

    index = _load_index(
        str(kb_root),
        max(1, int(chunk_size)),
        max(0, int(chunk_overlap)),
        tuple(exclude_list),
    )
    chunks = list(index.get("chunks") or [])
    if not chunks:
        return {"query": query, "results": []}

    q_tokens = _tokenize(query)
    scored = [(chunk, _score_chunk(q_tokens, chunk)) for chunk in chunks]
    scored.sort(key=lambda item: item[1], reverse=True)
    results: List[Dict[str, object]] = []
    for chunk, score in scored[: max(1, int(top_k))]:
        if score <= 0:
            continue
        entry = {
            "id": chunk.get("id"),
            "path": chunk.get("path"),
            "heading": chunk.get("heading"),
            "section_heading": chunk.get("section_heading"),
            "score": round(float(score), 3),
            "tags": chunk.get("tags") or [],
            "meta": chunk.get("meta") or {},
        }
        if include_text:
            text = str(chunk.get("text") or "")
            if max_chars > 0 and len(text) > max_chars:
                entry["text"] = text[: max_chars].rstrip() + "..."
                entry["truncated"] = True
            else:
                entry["text"] = text
                entry["truncated"] = False
        results.append(entry)

    return {"query": query, "results": results}


__all__ = ["search_knowledge_base"]
