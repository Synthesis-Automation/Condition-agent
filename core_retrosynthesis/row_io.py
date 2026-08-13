"""Canonical JSONL row iteration for retrosynthesis datasets."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, TextIO


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _row_files(
    source: str | Path,
    include: Iterable[str] = (),
) -> tuple[Path, ...]:
    path = Path(source)
    if path.is_file():
        return (path,)
    if not path.is_dir():
        raise FileNotFoundError(path)
    patterns = tuple(include) or ("*.jsonl", "*.jsonl.gz")
    return tuple(
        sorted(
            {
                candidate
                for pattern in patterns
                for candidate in path.rglob(pattern)
            },
            key=lambda value: value.as_posix(),
        )
    )


def iter_rows(
    source: str | Path,
    *,
    include: Iterable[str] = (),
) -> Iterator[Dict[str, Any]]:
    """Yield JSON object rows from one file or a directory tree."""

    for path in _row_files(source, include):
        with _open_text(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise ValueError(
                        f"{path}:{line_number}: row is not an object"
                    )
                yield value


__all__ = ["iter_rows"]
