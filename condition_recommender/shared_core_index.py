"""Atomic, source-bound projection artifacts for experimental shared retrieval.

This is a derived index over admitted canonical observations, not a separate
conversion path. Rebuilding it reuses stored graph evidence and never invents
unobserved reactant combinations.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import tempfile
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import closing
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterator

from reactive_taxonomy.shared_reaction_core import (
    SharedReactionCore,
    build_shared_reaction_core,
    shared_core_definition_hash,
)

from .generic_indexing import (
    GenericIndexedReaction,
    GenericReactionIndex,
    _indexed_reaction_payload,
)

STORAGE_VERSION = "1.0"


class SharedCoreBuildCancelled(RuntimeError):
    """A cancelled projection build leaves the previous artifact intact."""


def _encode(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _hash(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _row_digest(row: GenericIndexedReaction) -> str:
    return _hash(_encode(_indexed_reaction_payload(row)))


def _index_identity(index: GenericReactionIndex) -> str:
    identity = getattr(index.rows, "artifact_identity", None)
    if identity:
        return str(identity)
    digest = hashlib.sha256()
    for row in index.rows:
        digest.update(_row_digest(row).encode())
    return "SCBASE1:" + digest.hexdigest()


@dataclass(frozen=True)
class SharedCoreIndex:
    """Read-only projection artifact bound to one canonical precedent index."""

    path: Path
    source_identity: str
    definition_hash: str
    row_count: int

    def validate_binding(self, index: GenericReactionIndex) -> None:
        """Refuse mismatched observations, chemistry definitions or row order."""
        if (
            self.source_identity != _index_identity(index)
            or self.row_count != len(index.rows)
            or self.definition_hash != shared_core_definition_hash()
        ):
            raise ValueError(
                "SHARED_CORE_ARTIFACT_MISMATCH: rebuild for this precedent index"
            )

    def lookup(self, kind: str, key: str, limit: int) -> tuple[tuple[int, ...], bool]:
        """Return deterministic bounded positions and explicit truncation."""
        if limit < 1:
            raise ValueError("candidate limit must be positive")
        with closing(
            sqlite3.connect(f"{self.path.resolve().as_uri()}?mode=ro", uri=True)
        ) as connection:
            rows = connection.execute(
                "SELECT position FROM lookup WHERE kind=? AND key=? ORDER BY position LIMIT ?",
                (kind, key, limit + 1),
            ).fetchall()
        return tuple(row[0] for row in rows[:limit]), len(rows) > limit

    def projection(
        self, position: int, row: GenericIndexedReaction
    ) -> SharedReactionCore:
        """Validate precedent binding and payload integrity before comparison."""
        with closing(
            sqlite3.connect(f"{self.path.resolve().as_uri()}?mode=ro", uri=True)
        ) as connection:
            value = connection.execute(
                "SELECT row_hash, payload, payload_hash FROM projection WHERE position=?",
                (position,),
            ).fetchone()
        if value is None or value[0] != _row_digest(row) or value[2] != _hash(value[1]):
            raise ValueError(
                "SHARED_CORE_OBSERVATION_MISMATCH: rebuild projection artifact"
            )
        return SharedReactionCore.from_dict(json.loads(value[1]))


def load_shared_core_index(
    path: str | Path, index: GenericReactionIndex
) -> SharedCoreIndex:
    """Load an explicitly selected artifact; never silently use stale keys."""
    source = Path(path).resolve()
    with closing(sqlite3.connect(f"{source.as_uri()}?mode=ro", uri=True)) as connection:
        metadata = json.loads(
            connection.execute("SELECT payload FROM metadata").fetchone()[0]
        )
        count = connection.execute("SELECT count(*) FROM projection").fetchone()[0]
    if (
        metadata.get("storage_version") != STORAGE_VERSION
        or count != metadata.get("row_count")
        or metadata.get("build_complete", True) is not True
    ):
        raise ValueError("invalid shared core artifact manifest")
    result = SharedCoreIndex(
        source, metadata["source_identity"], metadata["definition_hash"], count
    )
    result.validate_binding(index)
    return result


def _project_rows(
    start: int,
    rows: tuple[GenericIndexedReaction, ...],
) -> list[
    tuple[
        int,
        str,
        str,
        str,
        tuple[tuple[str, str], ...],
        tuple[str, ...],
        tuple[str, ...],
    ]
]:
    """Calculate one bounded batch without sharing RDKit or SQLite state."""
    results = []
    for position, row in enumerate(rows, start):
        projection = build_shared_reaction_core(
            row.reaction_smiles, row.signature, row.reaction_core
        )
        payload = _encode(projection.to_dict())
        keys = set()
        if projection.levels:
            keys.update(
                (
                    ("whole_reaction", projection.reaction_identity),
                    ("product_identity", projection.product_identity),
                )
            )
        for level in projection.levels:
            keys.add((level.level, level.key))
            keys.add(("product_core", level.product_side_key))
        results.append(
            (
                position,
                _row_digest(row),
                payload,
                _hash(payload),
                tuple(sorted(keys)),
                tuple(level.level for level in projection.levels),
                projection.unavailable_reasons,
            )
        )
    return results


def _project_source_batch(source: str, total: int, start: int, stop: int) -> list:
    """Read and project a source range inside a worker, avoiding large IPC inputs."""
    from .sqlite_indexing import SQLiteReactionRows

    return _project_rows(
        start, SQLiteReactionRows(Path(source), total).select(range(start, stop))
    )


def _projection_batches(
    index: GenericReactionIndex,
    start: int,
    workers: int,
    batch_size: int,
    cancel_check: Callable[[], bool] | None,
) -> Iterator[list]:
    ranges = iter(
        (offset, min(offset + batch_size, len(index.rows)))
        for offset in range(start, len(index.rows), batch_size)
    )
    if workers == 1:
        for first, stop in ranges:
            if cancel_check is not None and cancel_check():
                raise SharedCoreBuildCancelled("Shared-core projection build cancelled")
            yield _project_rows(first, index.select(range(first, stop)))
        return
    source = str(index.rows.source_path)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pending = deque()

        def submit() -> None:
            item = next(ranges, None)
            if item is not None:
                pending.append(
                    executor.submit(
                        _project_source_batch, source, len(index.rows), *item
                    )
                )

        for _ in range(workers * 2):
            submit()
        while pending:
            if cancel_check is not None and cancel_check():
                for future in pending:
                    future.cancel()
                raise SharedCoreBuildCancelled("Shared-core projection build cancelled")
            yield pending.popleft().result()
            submit()


def build_shared_core_index(
    index: GenericReactionIndex,
    path: str | Path,
    *,
    cancel_check: Callable[[], bool] | None = None,
    progress_callback: Callable[[int], None] | None = None,
    workers: int = 1,
    batch_size: int = 128,
    resume: bool = False,
) -> dict[str, Any]:
    """Backfill atomically, optionally retaining source-bound batch checkpoints.

    Parallel workers read bounded source ranges; one writer commits in source
    order. Chemistry keys and row hashes are identical to a serial build.
    Resume checkpoints are never served as a completed projection artifact.
    """
    if type(workers) is not int or not 1 <= workers <= 32:
        raise ValueError("workers must be between 1 and 32")
    if type(batch_size) is not int or not 1 <= batch_size <= 1024:
        raise ValueError("batch_size must be between 1 and 1024")
    target = Path(path).resolve()
    source_path = getattr(index.rows, "source_path", None)
    if workers > 1 and source_path is None:
        raise ValueError("parallel projection builds require a persisted SQLite index")
    if source_path is not None and target == Path(source_path).resolve():
        raise ValueError("projection artifact must not overwrite its source index")
    target.parent.mkdir(parents=True, exist_ok=True)
    if resume:
        temporary = target.with_name(target.name + ".building")
    else:
        descriptor, name = tempfile.mkstemp(
            prefix=target.name + ".", suffix=".tmp", dir=target.parent
        )
        os.close(descriptor)
        temporary = Path(name)
    counts: Counter[str] = Counter()
    reasons: Counter[str] = Counter()
    metadata = dict(
        storage_version=STORAGE_VERSION,
        source_identity=_index_identity(index),
        definition_hash=shared_core_definition_hash(),
        row_count=len(index.rows),
        status="experimental_pending_independent_review",
        build_complete=False,
        completed_rows=0,
    )
    try:
        with closing(sqlite3.connect(temporary)) as connection, connection:
            connection.execute("PRAGMA cache_size=-32768")
            existing = connection.execute(
                "SELECT name FROM sqlite_master WHERE name='metadata'"
            ).fetchone()
            if existing:
                checkpoint = json.loads(
                    connection.execute("SELECT payload FROM metadata").fetchone()[0]
                )
                for key in (
                    "storage_version",
                    "source_identity",
                    "definition_hash",
                    "row_count",
                ):
                    if checkpoint.get(key) != metadata[key]:
                        raise ValueError(
                            "SHARED_CORE_CHECKPOINT_MISMATCH: use a matching source or new output path"
                        )
                completed = int(checkpoint.get("completed_rows", 0))
                count, last = connection.execute(
                    "SELECT count(*), max(position) FROM projection"
                ).fetchone()
                if count != completed or (completed and last != completed - 1):
                    raise ValueError("invalid shared core checkpoint prefix")
                metadata.update(checkpoint)
                counts.update(checkpoint.get("eligible_counts", {}))
                reasons.update(checkpoint.get("unavailable_reasons", {}))
            else:
                connection.executescript("""
                CREATE TABLE metadata (payload TEXT NOT NULL);
                CREATE TABLE projection (position INTEGER PRIMARY KEY, row_hash TEXT NOT NULL,
                    payload TEXT NOT NULL, payload_hash TEXT NOT NULL);
                CREATE TABLE lookup (kind TEXT NOT NULL, key TEXT NOT NULL, position INTEGER NOT NULL,
                    PRIMARY KEY (kind, key, position), FOREIGN KEY(position) REFERENCES projection(position))
                    WITHOUT ROWID;
                """)
                connection.execute(
                    "INSERT INTO metadata VALUES (?)", (_encode(metadata),)
                )
                connection.commit()
            if progress_callback is not None:
                progress_callback(metadata["completed_rows"])
            for batch in _projection_batches(
                index, metadata["completed_rows"], workers, batch_size, cancel_check
            ):
                if cancel_check is not None and cancel_check():
                    raise SharedCoreBuildCancelled(
                        "Shared-core projection build cancelled"
                    )
                for (
                    position,
                    row_hash,
                    payload,
                    payload_hash,
                    keys,
                    levels,
                    unavailable,
                ) in batch:
                    connection.execute(
                        "INSERT INTO projection VALUES (?, ?, ?, ?)",
                        (position, row_hash, payload, payload_hash),
                    )
                    connection.executemany(
                        "INSERT INTO lookup VALUES (?, ?, ?)",
                        [(kind, key, position) for kind, key in keys],
                    )
                    counts.update(levels)
                    reasons.update(unavailable)
                metadata.update(
                    completed_rows=batch[-1][0] + 1,
                    eligible_counts=dict(sorted(counts.items())),
                    unavailable_reasons=dict(sorted(reasons.items())),
                )
                connection.execute(
                    "UPDATE metadata SET payload=?", (_encode(metadata),)
                )
                connection.commit()
                if progress_callback is not None:
                    progress_callback(metadata["completed_rows"])
            metadata.update(
                eligible_counts=dict(sorted(counts.items())),
                unavailable_reasons=dict(sorted(reasons.items())),
            )
            if _index_identity(index) != metadata["source_identity"]:
                raise ValueError("source index changed during projection build")
            if connection.execute("PRAGMA foreign_key_check").fetchall():
                raise ValueError("dangling shared core lookup")
            metadata["build_complete"] = True
            connection.execute("UPDATE metadata SET payload=?", (_encode(metadata),))
        os.replace(temporary, target)
    finally:
        if not resume:
            temporary.unlink(missing_ok=True)
    if progress_callback is not None:
        progress_callback(len(index.rows))
    return metadata
