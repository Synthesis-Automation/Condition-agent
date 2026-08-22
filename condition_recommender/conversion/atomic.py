"""Atomic filesystem writes shared by conversion artifact builders."""

from __future__ import annotations

import json
import os
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Mapping

_REPLACE_RETRY_DELAYS_SECONDS = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6)


def replace_with_retry(source: Path, destination: Path) -> None:
    """Atomically replace ``destination``, tolerating brief Windows locks.

    Virus scanners, indexers, and readers can transiently prevent a replace on
    Windows. Retrying ``PermissionError`` is safe because the source remains in
    place after a failed replacement. Other filesystem errors are not masked.
    """
    for delay in (*_REPLACE_RETRY_DELAYS_SECONDS, None):
        try:
            os.replace(source, destination)
            return
        except PermissionError:
            if delay is None:
                raise
            time.sleep(delay)


@contextmanager
def atomic_output_path(destination: str | Path) -> Iterator[Path]:
    """Yield a unique sibling path and atomically publish it on success."""
    output = Path(destination)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        yield temporary
        replace_with_retry(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_json(path: str | Path, payload: Mapping[str, Any]) -> None:
    """Write a deterministic JSON object without exposing a partial file."""
    with atomic_output_path(path) as temporary:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
