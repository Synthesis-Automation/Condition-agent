"""
Extract condition summary tables from the curated knowledge base.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Sequence


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

_TABLE_SEPARATOR_RE = re.compile(r"^[\s|\-:]+$")
_COLUMN_MAP = {
    "catalyst": "catalyst",
    "ligand": "ligand",
    "base": "base",
    "solvent": "solvent",
    "secondary solvent": "secondary_solvent",
    "additive": "additive",
    "coupling reagent": "coupling_reagent",
    "temp (c)": "temperature_c",
    "temp": "temperature_c",
    "temperature": "temperature_c",
    "time (h)": "time_h",
    "time": "time_h",
    "notes": "notes",
    "rank": "rank",
    "reaction id": "reaction_id",
    "reaction type": "reaction_type",
}


def _tokenize(text: str) -> List[str]:
    tokens = re.findall(r"[A-Za-z0-9\\-]+", text.lower())
    return [t for t in tokens if t and t not in _STOPWORDS]


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="utf-8", errors="replace")


def _parse_front_matter(lines: Sequence[str]) -> tuple[Dict[str, str], int]:
    if not lines or lines[0].strip() != "---":
        return {}, 0
    end_idx = None
    for idx in range(1, len(lines)):
        if lines[idx].strip() == "---":
            end_idx = idx
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


def _parse_tags(meta: Dict[str, str]) -> List[str]:
    raw = meta.get("tags", "")
    tags = [tag.strip() for tag in raw.split(",") if tag.strip()]
    return tags


def _normalize_column_name(name: str) -> str:
    cleaned = re.sub(r"\s+", " ", name.strip().lower())
    return cleaned


def _split_table_row(line: str) -> List[str]:
    parts = [part.strip() for part in line.strip().strip("|").split("|")]
    return parts


def _iter_markdown_tables(lines: Sequence[str]) -> Iterable[Dict[str, object]]:
    current_heading = ""
    current_label = ""
    idx = 0
    while idx < len(lines):
        line = lines[idx].rstrip()
        stripped = line.strip()

        if stripped.startswith("#"):
            current_heading = stripped.lstrip("#").strip()
            current_label = ""
            idx += 1
            continue

        if stripped.endswith(":") and stripped and not stripped.startswith("|"):
            current_label = stripped.rstrip(":").strip()
            idx += 1
            continue

        if stripped.startswith("|") and idx + 1 < len(lines):
            header_line = lines[idx].strip()
            separator_line = lines[idx + 1].strip()
            if not separator_line.startswith("|") or not _TABLE_SEPARATOR_RE.match(separator_line):
                idx += 1
                continue
            headers = _split_table_row(header_line)
            idx += 2
            rows: List[List[str]] = []
            while idx < len(lines):
                row_line = lines[idx].strip()
                if not row_line or not row_line.startswith("|"):
                    break
                rows.append(_split_table_row(row_line))
                idx += 1
            context_parts = [part for part in (current_heading, current_label) if part]
            yield {
                "headers": headers,
                "rows": rows,
                "heading": current_heading,
                "label": current_label,
                "context": " | ".join(context_parts),
            }
            continue

        idx += 1


@dataclass
class _ConditionRow:
    condition: Dict[str, str]
    extras: Dict[str, str]
    context: str
    heading: str
    label: str


def _build_condition_rows(table: Dict[str, object]) -> List[_ConditionRow]:
    headers = [str(h) for h in table.get("headers") or []]
    rows = table.get("rows") or []
    context = str(table.get("context") or "")
    heading = str(table.get("heading") or "")
    label = str(table.get("label") or "")
    output: List[_ConditionRow] = []
    for raw_row in rows:
        row_values = [str(cell) for cell in raw_row]
        if len(row_values) < len(headers):
            row_values.extend([""] * (len(headers) - len(row_values)))
        row_map = dict(zip(headers, row_values))
        condition: Dict[str, str] = {}
        extras: Dict[str, str] = {}
        for header, value in row_map.items():
            norm = _normalize_column_name(header)
            key = _COLUMN_MAP.get(norm)
            if key:
                condition[key] = value
            else:
                extras[header] = value
        if not condition:
            continue
        output.append(
            _ConditionRow(
                condition=condition,
                extras=extras,
                context=context,
                heading=heading,
                label=label,
            )
        )
    return output


def _iter_markdown_files(root: Path) -> Iterable[Path]:
    if root.is_file():
        yield root
        return
    for path in root.rglob("*.md"):
        if path.is_file():
            yield path


@lru_cache(maxsize=4)
def _load_condition_rows(root: str) -> List[Dict[str, object]]:
    root_path = Path(root)
    if not root_path.exists():
        return []

    entries: List[Dict[str, object]] = []
    for path in sorted(_iter_markdown_files(root_path)):
        text = _read_text(path)
        lines = text.splitlines()
        meta, body_start = _parse_front_matter(lines)
        tags = _parse_tags(meta)
        title = meta.get("title", "")
        doc_id = meta.get("id", "")
        scope = meta.get("scope", "")
        body_lines = lines[body_start:] if body_start else lines
        for table in _iter_markdown_tables(body_lines):
            rows = _build_condition_rows(table)
            if not rows:
                continue
            context = " | ".join(
                part
                for part in (
                    title,
                    doc_id,
                    scope,
                    table.get("context"),
                )
                if part
            )
            for row in rows:
                entries.append(
                    {
                        "condition": row.condition,
                        "extras": row.extras,
                        "context": context,
                        "heading": row.heading,
                        "label": row.label,
                        "tags": tags,
                        "doc_id": doc_id,
                        "path": path.as_posix(),
                    }
                )
    return entries


def _score_entry(query_tokens: List[str], entry: Dict[str, object]) -> float:
    context = str(entry.get("context") or "")
    tokens = _tokenize(context)
    overlap = len(set(query_tokens) & set(tokens))
    tag_bonus = 0.0
    for tag in entry.get("tags") or []:
        if str(tag).lower() in query_tokens:
            tag_bonus += 1.5
    doc_id = str(entry.get("doc_id") or "").lower()
    if doc_id and doc_id in query_tokens:
        tag_bonus += 2.0
    return overlap + tag_bonus


def search_condition_summaries(
    query: str,
    *,
    top_k: int = 5,
    root: Optional[str] = None,
) -> Dict[str, object]:
    if not query or not query.strip():
        return {"query": query, "results": []}

    repo_root = Path(__file__).resolve().parents[1]
    kb_root = Path(root) if root else (repo_root / "data" / "knowledge_base")
    entries = _load_condition_rows(str(kb_root))
    if not entries:
        return {"query": query, "results": []}

    q_tokens = _tokenize(query)
    scored: List[Dict[str, object]] = []
    for entry in entries:
        score = _score_entry(q_tokens, entry)
        if score <= 0:
            continue
        scored.append(
            {
                "score": round(float(score), 3),
                "condition": entry.get("condition") or {},
                "extras": entry.get("extras") or {},
                "context": entry.get("context") or "",
                "heading": entry.get("heading") or "",
                "label": entry.get("label") or "",
                "tags": entry.get("tags") or [],
                "doc_id": entry.get("doc_id") or "",
                "path": entry.get("path") or "",
            }
        )

    scored.sort(key=lambda item: item["score"], reverse=True)
    results = []
    for idx, item in enumerate(scored[: max(1, int(top_k))], start=1):
        payload = dict(item)
        payload["rank"] = idx
        results.append(payload)

    return {"query": query, "results": results}


__all__ = ["search_condition_summaries"]
