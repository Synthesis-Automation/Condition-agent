"""Lazy SQLite persistence for the interactive recommendation index."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import zlib
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Tuple, overload

from reactive_taxonomy import departing_fragment_tokens

from .edit_prototypes import AnonymousEditPrototype, anonymous_edit_prototype
from .fallback_similarity import fallback_index_tokens
from .reaction_facets import reaction_facet_keys
from .generic_indexing import (
    GENERIC_INDEX_SCHEMA_VERSION,
    GenericIndexedReaction,
    GenericReactionIndex,
    _iter_generic_index_rows,
    _index_maps,
    _indexed_reaction_from_payload,
    _indexed_reaction_payload,
    _validate_index_rows,
    _validate_index_metadata,
    build_generic_index_from_rows,
)
from .models import (
    FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION,
    PrecedentIndexScope,
    PrecedentTier,
)
from .signature_features import environment_tokens


SQLITE_INDEX_STORAGE_SCHEMA_VERSION = "1.0"
_ROW_ENCODING = "json+zlib.v1"
_POSITIONS_ENCODING = "json+zlib.v1"


def _indexed_edit_prototype(row: GenericIndexedReaction) -> AnonymousEditPrototype | None:
    """Project cycles wholly contained in departing fragments before indexing."""

    prototype = anonymous_edit_prototype(row.signature)
    if prototype is None or prototype.ring_count_delta >= 0:
        return prototype
    ring_count = 0
    for token in departing_fragment_tokens(row.reaction_smiles, row.signature):
        parts = token.split(":", 3)
        if len(parts) >= 3 and parts[1].startswith("R"):
            ring_count += int(parts[1][1:] or 0)
    if not ring_count:
        return prototype
    return replace(
        prototype,
        ring_count_delta=min(0, prototype.ring_count_delta + ring_count),
    )


class SQLiteIndexBuildCancelled(RuntimeError):
    """Raised after safely abandoning an incomplete streamed SQLite index."""


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _compress_json(value: Any) -> bytes:
    return zlib.compress(_canonical_json_bytes(value), level=6)


def _decompress_json(value: bytes) -> Any:
    return json.loads(zlib.decompress(value))


def _read_connection(path: Path) -> sqlite3.Connection:
    uri = path.resolve().as_uri() + "?mode=ro"
    return sqlite3.connect(uri, uri=True, timeout=30.0)


class SQLiteReactionRows(Sequence[GenericIndexedReaction]):
    """Read-only sequence that materializes precedent rows on demand."""

    __hash__ = object.__hash__

    def __init__(self, path: Path, row_count: int) -> None:
        self._path = path
        self._row_count = row_count

    def __len__(self) -> int:
        return self._row_count

    @overload
    def __getitem__(self, position: int) -> GenericIndexedReaction: ...

    @overload
    def __getitem__(self, position: slice) -> Sequence[GenericIndexedReaction]: ...

    def __getitem__(
        self,
        position: int | slice,
    ) -> GenericIndexedReaction | Sequence[GenericIndexedReaction]:
        if isinstance(position, slice):
            return self.select(range(*position.indices(self._row_count)))
        normalized = position + self._row_count if position < 0 else position
        if normalized < 0 or normalized >= self._row_count:
            raise IndexError(position)
        return self._row(normalized)

    @lru_cache(maxsize=4_096)
    def _row(self, position: int) -> GenericIndexedReaction:
        with _read_connection(self._path) as connection:
            result = connection.execute(
                "SELECT payload FROM precedents WHERE position = ?",
                (position,),
            ).fetchone()
        if result is None:
            raise ValueError(f"SQLite index is missing precedent row {position}")
        return _indexed_reaction_from_payload(_decompress_json(result[0]))

    def __iter__(self) -> Iterator[GenericIndexedReaction]:
        with _read_connection(self._path) as connection:
            cursor = connection.execute(
                "SELECT position, payload FROM precedents ORDER BY position"
            )
            expected = 0
            for position, payload in cursor:
                if int(position) != expected:
                    raise ValueError(
                        "SQLite index precedent positions are not contiguous"
                    )
                expected += 1
                yield _indexed_reaction_from_payload(_decompress_json(payload))
        if expected != self._row_count:
            raise ValueError("SQLite index precedent row count mismatch")

    def select(
        self,
        positions: Iterable[int],
    ) -> Tuple[GenericIndexedReaction, ...]:
        """Fetch an ordered candidate subset with bounded SQL statements."""
        requested = tuple(int(value) for value in positions)
        if not requested:
            return ()
        for position in requested:
            if position < 0 or position >= self._row_count:
                raise IndexError(position)
        loaded: Dict[int, GenericIndexedReaction] = {}
        unique_positions = tuple(sorted(set(requested)))
        with _read_connection(self._path) as connection:
            for start in range(0, len(unique_positions), 500):
                chunk = unique_positions[start : start + 500]
                placeholders = ",".join("?" for _ in chunk)
                rows = connection.execute(
                    "SELECT position, payload FROM precedents "
                    f"WHERE position IN ({placeholders})",
                    chunk,
                )
                for position, payload in rows:
                    loaded[int(position)] = _indexed_reaction_from_payload(
                        _decompress_json(payload)
                    )
        missing = set(requested) - set(loaded)
        if missing:
            raise ValueError(
                "SQLite index is missing precedent row(s): "
                + ", ".join(str(value) for value in sorted(missing))
            )
        return tuple(loaded[position] for position in requested)

    def edit_graph_candidate_positions(
        self,
        prototype: AnonymousEditPrototype,
    ) -> Tuple[int, ...]:
        """Return the SQLite-prefiltered superset for edit-graph scoring."""
        ring_sign = (prototype.ring_count_delta > 0) - (
            prototype.ring_count_delta < 0
        )
        keys = tuple(
            f"{ring_sign}:{pair}"
            for pair in sorted(set(prototype.formed_element_pairs))
        )
        if not keys:
            return ()
        placeholders = ",".join("?" for _ in keys)
        positions = set()
        with _read_connection(self._path) as connection:
            rows = connection.execute(
                "SELECT positions FROM lookups WHERE lookup_kind = ? "
                f"AND lookup_key IN ({placeholders})",
                ("edit_graph_features", *keys),
            )
            for (payload,) in rows:
                positions.update(int(value) for value in _decompress_json(payload))
        return tuple(sorted(positions))


class SQLiteLookupMapping(Mapping[str, Tuple[int, ...]]):
    """Read-only lookup map backed by one SQLite lookup category."""

    __hash__ = object.__hash__

    def __init__(self, path: Path, kind: str, key_count: int) -> None:
        self._path = path
        self._kind = kind
        self._key_count = key_count

    def __len__(self) -> int:
        return self._key_count

    def __iter__(self) -> Iterator[str]:
        with _read_connection(self._path) as connection:
            rows = connection.execute(
                "SELECT lookup_key FROM lookups "
                "WHERE lookup_kind = ? ORDER BY lookup_key",
                (self._kind,),
            )
            yield from (str(row[0]) for row in rows)

    @lru_cache(maxsize=4_096)
    def __getitem__(self, key: str) -> Tuple[int, ...]:
        normalized = str(key)
        with _read_connection(self._path) as connection:
            result = connection.execute(
                "SELECT positions FROM lookups "
                "WHERE lookup_kind = ? AND lookup_key = ?",
                (self._kind, normalized),
            ).fetchone()
        if result is None:
            raise KeyError(normalized)
        return tuple(int(value) for value in _decompress_json(result[0]))


def _metadata_payload(
    index: GenericReactionIndex,
    *,
    lookup_counts: Mapping[str, int],
    auxiliary_lookup_counts: Mapping[str, int],
    index_id: str,
) -> Dict[str, Any]:
    return {
        "artifact_type": "generic_reaction_index",
        "schema_version": GENERIC_INDEX_SCHEMA_VERSION,
        "storage_schema_version": SQLITE_INDEX_STORAGE_SCHEMA_VERSION,
        "row_encoding": _ROW_ENCODING,
        "positions_encoding": _POSITIONS_ENCODING,
        "reaction_signature_schema_version": index.reaction_signature_schema_version,
        "reaction_core_schema_version": index.reaction_core_schema_version,
        "reaction_core_algorithm_version": index.reaction_core_algorithm_version,
        "taxonomy_definition_versions": dict(index.taxonomy_definition_versions),
        "fallback_descriptor_schema_version": (
            index.fallback_descriptor_schema_version
        ),
        "fallback_definition_versions": dict(index.fallback_definition_versions),
        "record_schema_versions": list(index.record_schema_versions),
        "converter_definition_versions": list(index.converter_definition_versions),
        "precedent_scope": index.precedent_scope.value,
        "core_eligibility_definition_version": (
            index.core_eligibility_definition_version
        ),
        "fragment_source_capability_definition_version": (
            FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION
        ),
        "index_id": index_id,
        "row_count": len(index.rows),
        "lookup_counts": dict(sorted(lookup_counts.items())),
        "auxiliary_lookup_counts": dict(sorted(auxiliary_lookup_counts.items())),
    }


def build_sqlite_generic_index(
    records: Iterable[Mapping[str, Any]],
    path: str | Path,
    *,
    include_review: bool = False,
    progress_callback: Callable[[str, int], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> Dict[str, Any]:
    """Build a deterministic SQLite index with bounded Python memory.

    Admitted rows and lookup memberships are staged in a temporary SQLite
    database. Sorting and aggregation therefore spill to disk instead of
    retaining the full corpus, row payloads, and lookup maps in Python.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    staging = destination.with_suffix(destination.suffix + ".stage")
    for transient in (temporary, staging):
        if transient.is_file():
            transient.unlink()

    scope = (
        PrecedentIndexScope.TRUSTED_AND_REVIEW_CORE
        if include_review
        else PrecedentIndexScope.TRUSTED
    )
    stage_connection: sqlite3.Connection | None = None
    final_connection: sqlite3.Connection | None = None
    completed = False
    try:
        stage_connection = sqlite3.connect(staging)
        stage_connection.executescript(
            """
            PRAGMA journal_mode = OFF;
            PRAGMA synchronous = OFF;
            PRAGMA temp_store = FILE;
            PRAGMA cache_size = -32768;
            PRAGMA page_size = 32768;
            CREATE TABLE staged_rows (
                ordinal INTEGER PRIMARY KEY,
                canonical_reaction_id TEXT NOT NULL,
                reaction_id TEXT NOT NULL,
                observation_id TEXT NOT NULL,
                recipe_id TEXT NOT NULL,
                payload BLOB NOT NULL
            );
            CREATE TABLE lookup_entries (
                lookup_kind TEXT NOT NULL,
                lookup_key TEXT NOT NULL,
                position INTEGER NOT NULL,
                PRIMARY KEY (lookup_kind, lookup_key, position)
            ) WITHOUT ROWID;
            """
        )
        stage_batch = []
        eligible_count = 0
        for row in _iter_generic_index_rows(records, include_review=include_review):
            if cancel_check is not None and cancel_check():
                raise SQLiteIndexBuildCancelled(
                    "Index build cancelled while scanning canonical records"
                )
            _validate_index_rows((row,))
            if scope == PrecedentIndexScope.TRUSTED:
                if row.precedent_tier != PrecedentTier.TRUSTED:
                    raise ValueError(
                        "Trusted index cannot contain review-core precedents"
                    )
            elif row.precedent_tier not in {
                PrecedentTier.TRUSTED,
                PrecedentTier.REVIEW_CORE,
            }:
                raise ValueError("Unsupported precedent tier")
            canonical = _canonical_json_bytes(_indexed_reaction_payload(row))
            stage_batch.append(
                (
                    eligible_count,
                    row.canonical_reaction_id,
                    row.reaction_id,
                    row.observation_id,
                    row.recipe_id,
                    zlib.compress(canonical, level=6),
                )
            )
            eligible_count += 1
            if len(stage_batch) >= 100:
                stage_connection.executemany(
                    "INSERT INTO staged_rows VALUES (?, ?, ?, ?, ?, ?)",
                    stage_batch,
                )
                stage_batch.clear()
            if eligible_count % 1_000 == 0 and progress_callback is not None:
                progress_callback("staged", eligible_count)
        if stage_batch:
            stage_connection.executemany(
                "INSERT INTO staged_rows VALUES (?, ?, ?, ?, ?, ?)",
                stage_batch,
            )
        stage_connection.commit()

        final_connection = sqlite3.connect(temporary)
        final_connection.executescript(
            """
            PRAGMA journal_mode = OFF;
            PRAGMA synchronous = OFF;
            PRAGMA temp_store = FILE;
            PRAGMA cache_size = -32768;
            PRAGMA page_size = 32768;
            CREATE TABLE metadata (
                singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
                payload TEXT NOT NULL
            );
            CREATE TABLE precedents (
                position INTEGER PRIMARY KEY,
                payload BLOB NOT NULL
            );
            CREATE TABLE lookups (
                lookup_kind TEXT NOT NULL,
                lookup_key TEXT NOT NULL,
                positions BLOB NOT NULL,
                PRIMARY KEY (lookup_kind, lookup_key)
            ) WITHOUT ROWID;
            """
        )
        digest = hashlib.sha256()
        row_batch = []
        lookup_batch = []
        sorted_rows = stage_connection.execute(
            "SELECT payload FROM staged_rows ORDER BY "
            "canonical_reaction_id, reaction_id, observation_id, recipe_id, ordinal"
        )
        for position, (compressed_payload,) in enumerate(sorted_rows):
            if cancel_check is not None and cancel_check():
                raise SQLiteIndexBuildCancelled(
                    "Index build cancelled while writing precedent rows"
                )
            payload = _decompress_json(compressed_payload)
            row = _indexed_reaction_from_payload(payload)
            canonical = _canonical_json_bytes(payload)
            digest.update(b"row\0")
            digest.update(str(position).encode("ascii"))
            digest.update(b"\0")
            digest.update(canonical)
            row_batch.append((position, compressed_payload))

            signature_fields = {
                "exact": "exact_signature_key",
                "handles": "handle_signature_key",
                "transformations": "transformation_signature_key",
                "bond_edits": "bond_edit_signature_key",
                "environments": "environment_signature_key",
            }
            for kind, field in signature_fields.items():
                key = str(row.signature.get(field) or "")
                if key:
                    lookup_batch.append((kind, key, position))
            core_fields = {
                "core_exact": "exact_core_key",
                "core_typed": "typed_core_key",
                "core_shapes": "shape_core_key",
                "core_centers": "center_transition_key",
            }
            for kind, field in core_fields.items():
                key = str(row.reaction_core.get(field) or "")
                if key:
                    lookup_batch.append((kind, key, position))
            facet_keys = reaction_facet_keys(
                row.signature,
                row.reaction_core,
                row.fallback_descriptor,
            )
            exact_facet = facet_keys.get("reaction_facet_exact")
            if exact_facet:
                lookup_batch.append(("facet_exact", exact_facet, position))
            attachment_facet = facet_keys.get(
                "reaction_facet_attachment_relaxed"
            )
            if attachment_facet:
                lookup_batch.append(
                    ("facet_attachments", attachment_facet, position)
                )
            if row.named_family:
                lookup_batch.append(("families", row.named_family, position))
            for token in set(environment_tokens(row.signature)):
                lookup_batch.append(("environment_features", token, position))
            fragment_tokens = set(
                departing_fragment_tokens(row.reaction_smiles, row.signature)
            )
            if fragment_tokens:
                lookup_batch.append(("departing_fragments", "DF1:KNOWN", position))
            for token in fragment_tokens:
                lookup_batch.append(("departing_fragments", token, position))
            for token in fallback_index_tokens(row.fallback_descriptor):
                lookup_batch.append(("fallback_features", token, position))
            partial_key = str(
                row.fallback_descriptor.get("partial_transformation_key") or ""
            )
            if partial_key:
                lookup_batch.append(
                    ("partial_transformations", partial_key, position)
                )
            prototype = _indexed_edit_prototype(row)
            if prototype is not None:
                ring_sign = (prototype.ring_count_delta > 0) - (
                    prototype.ring_count_delta < 0
                )
                for pair in sorted(set(prototype.formed_element_pairs)):
                    lookup_batch.append(
                        ("edit_graph_features", f"{ring_sign}:{pair}", position)
                    )

            if len(row_batch) >= 100:
                final_connection.executemany(
                    "INSERT INTO precedents(position, payload) VALUES (?, ?)",
                    row_batch,
                )
                row_batch.clear()
            if len(lookup_batch) >= 5_000:
                stage_connection.executemany(
                    "INSERT OR IGNORE INTO lookup_entries VALUES (?, ?, ?)",
                    lookup_batch,
                )
                lookup_batch.clear()
            if (position + 1) % 1_000 == 0 and progress_callback is not None:
                progress_callback("rows", position + 1)
        if row_batch:
            final_connection.executemany(
                "INSERT INTO precedents(position, payload) VALUES (?, ?)",
                row_batch,
            )
        if lookup_batch:
            stage_connection.executemany(
                "INSERT OR IGNORE INTO lookup_entries VALUES (?, ?, ?)",
                lookup_batch,
            )
        stage_connection.commit()

        lookup_names = (
            "exact",
            "handles",
            "transformations",
            "bond_edits",
            "environments",
            "facet_exact",
            "facet_attachments",
            "core_exact",
            "core_typed",
            "core_shapes",
            "core_centers",
            "environment_features",
            "departing_fragments",
            "fallback_features",
            "partial_transformations",
            "families",
        )
        lookup_counts = {name: 0 for name in lookup_names}
        edit_graph_key_order: Dict[str, int] = {}
        current_kind = ""
        current_key = ""
        positions: list[int] = []
        written_lookup_count = 0
        final_lookup_batch = []

        def write_lookup() -> None:
            nonlocal written_lookup_count
            if not current_kind:
                return
            canonical_positions = _canonical_json_bytes(tuple(positions))
            final_lookup_batch.append(
                (
                    current_kind,
                    current_key,
                    zlib.compress(canonical_positions, level=6),
                )
            )
            if len(final_lookup_batch) >= 2_000:
                final_connection.executemany(
                    "INSERT INTO lookups"
                    "(lookup_kind, lookup_key, positions) VALUES (?, ?, ?)",
                    final_lookup_batch,
                )
                final_lookup_batch.clear()
            if current_kind == "edit_graph_features":
                edit_graph_key_order[current_key] = positions[0]
            else:
                lookup_counts[current_kind] += 1
            written_lookup_count += 1
            if (
                written_lookup_count % 1_000 == 0
                and progress_callback is not None
            ):
                progress_callback("lookups", written_lookup_count)

        for kind, key, position in stage_connection.execute(
            "SELECT lookup_kind, lookup_key, position FROM lookup_entries "
            "ORDER BY lookup_kind, lookup_key, position"
        ):
            if cancel_check is not None and cancel_check():
                raise SQLiteIndexBuildCancelled(
                    "Index build cancelled while aggregating lookup keys"
                )
            if (kind, key) != (current_kind, current_key):
                write_lookup()
                current_kind = str(kind)
                current_key = str(key)
                positions = []
            positions.append(int(position))
        write_lookup()
        if final_lookup_batch:
            final_connection.executemany(
                "INSERT INTO lookups"
                "(lookup_kind, lookup_key, positions) VALUES (?, ?, ?)",
                final_lookup_batch,
            )

        for kind in sorted((*lookup_names, "edit_graph_features")):
            lookup_rows = final_connection.execute(
                "SELECT lookup_key, positions FROM lookups "
                "WHERE lookup_kind = ? ORDER BY lookup_key",
                (kind,),
            )
            if kind == "edit_graph_features":
                ordered_lookup_rows = sorted(
                    lookup_rows,
                    key=lambda item: (
                        edit_graph_key_order[str(item[0])],
                        str(item[0]),
                    ),
                )
            else:
                ordered_lookup_rows = lookup_rows
            for key, compressed_positions in ordered_lookup_rows:
                canonical_positions = zlib.decompress(compressed_positions)
                digest.update(b"lookup\0")
                digest.update(kind.encode("utf-8"))
                digest.update(b"\0")
                digest.update(str(key).encode("utf-8"))
                digest.update(b"\0")
                digest.update(canonical_positions)

        template_index = build_generic_index_from_rows(
            (),
            precedent_scope=scope,
        )
        metadata_index = replace(template_index, rows=range(eligible_count))
        metadata_identity = _metadata_payload(
            metadata_index,
            lookup_counts=lookup_counts,
            auxiliary_lookup_counts={
                "edit_graph_features": len(edit_graph_key_order),
            },
            index_id="",
        )
        digest.update(b"metadata\0")
        digest.update(_canonical_json_bytes(metadata_identity))
        metadata = _metadata_payload(
            metadata_index,
            lookup_counts=lookup_counts,
            auxiliary_lookup_counts={
                "edit_graph_features": len(edit_graph_key_order),
            },
            index_id="GRIS1:" + digest.hexdigest(),
        )
        final_connection.execute(
            "INSERT INTO metadata(singleton, payload) VALUES (1, ?)",
            (_canonical_json_bytes(metadata).decode("utf-8"),),
        )
        final_connection.commit()
        final_connection.close()
        final_connection = None
        stage_connection.close()
        stage_connection = None
        temporary.replace(destination)
        completed = True
    finally:
        if final_connection is not None:
            final_connection.close()
        if stage_connection is not None:
            stage_connection.close()
        if staging.is_file():
            staging.unlink()
        if not completed and temporary.is_file():
            temporary.unlink()

    return {
        "schema_version": metadata["schema_version"],
        "storage_schema_version": metadata["storage_schema_version"],
        "precedent_scope": metadata["precedent_scope"],
        "index_id": metadata["index_id"],
        "row_count": metadata["row_count"],
        "lookup_counts": metadata["lookup_counts"],
        "path": str(destination),
    }


def save_sqlite_generic_index(
    index: GenericReactionIndex,
    path: str | Path,
    *,
    index_id: str | None = None,
) -> Dict[str, Any]:
    """Atomically write a compact SQLite index for lazy interactive access."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    if temporary.is_file():
        temporary.unlink()
    lookup_maps = _index_maps(index)
    lookup_counts = {name: len(mapping) for name, mapping in lookup_maps.items()}
    edit_graph_features: Dict[str, list[int]] = defaultdict(list)
    digest = hashlib.sha256()
    completed = False
    connection = sqlite3.connect(temporary)
    try:
        connection.executescript(
            """
            PRAGMA journal_mode = OFF;
            PRAGMA synchronous = OFF;
            PRAGMA temp_store = MEMORY;
            PRAGMA page_size = 32768;
            CREATE TABLE metadata (
                singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
                payload TEXT NOT NULL
            );
            CREATE TABLE precedents (
                position INTEGER PRIMARY KEY,
                payload BLOB NOT NULL
            );
            CREATE TABLE lookups (
                lookup_kind TEXT NOT NULL,
                lookup_key TEXT NOT NULL,
                positions BLOB NOT NULL,
                PRIMARY KEY (lookup_kind, lookup_key)
            ) WITHOUT ROWID;
            """
        )
        row_batch = []
        for position, row in enumerate(index.rows):
            payload = _indexed_reaction_payload(row)
            prototype = _indexed_edit_prototype(row)
            if prototype is not None:
                ring_sign = (prototype.ring_count_delta > 0) - (
                    prototype.ring_count_delta < 0
                )
                for pair in sorted(set(prototype.formed_element_pairs)):
                    edit_graph_features[f"{ring_sign}:{pair}"].append(position)
            canonical = _canonical_json_bytes(payload)
            digest.update(b"row\0")
            digest.update(str(position).encode("ascii"))
            digest.update(b"\0")
            digest.update(canonical)
            row_batch.append((position, zlib.compress(canonical, level=6)))
            if len(row_batch) >= 500:
                connection.executemany(
                    "INSERT INTO precedents(position, payload) VALUES (?, ?)",
                    row_batch,
                )
                row_batch.clear()
        if row_batch:
            connection.executemany(
                "INSERT INTO precedents(position, payload) VALUES (?, ?)",
                row_batch,
            )
        lookup_batch = []
        persisted_maps = {
            **lookup_maps,
            "edit_graph_features": edit_graph_features,
        }
        for kind, mapping in sorted(persisted_maps.items()):
            for key, positions in mapping.items():
                values = tuple(int(value) for value in positions)
                canonical = _canonical_json_bytes(values)
                digest.update(b"lookup\0")
                digest.update(kind.encode("utf-8"))
                digest.update(b"\0")
                digest.update(str(key).encode("utf-8"))
                digest.update(b"\0")
                digest.update(canonical)
                lookup_batch.append(
                    (kind, str(key), zlib.compress(canonical, level=6))
                )
                if len(lookup_batch) >= 2_000:
                    connection.executemany(
                        "INSERT INTO lookups"
                        "(lookup_kind, lookup_key, positions) VALUES (?, ?, ?)",
                        lookup_batch,
                    )
                    lookup_batch.clear()
        if lookup_batch:
            connection.executemany(
                "INSERT INTO lookups"
                "(lookup_kind, lookup_key, positions) VALUES (?, ?, ?)",
                lookup_batch,
            )
        metadata_identity = _metadata_payload(
            index,
            lookup_counts=lookup_counts,
            auxiliary_lookup_counts={
                "edit_graph_features": len(edit_graph_features),
            },
            index_id="",
        )
        digest.update(b"metadata\0")
        digest.update(_canonical_json_bytes(metadata_identity))
        resolved_index_id = index_id or "GRIS1:" + digest.hexdigest()
        metadata = _metadata_payload(
            index,
            lookup_counts=lookup_counts,
            auxiliary_lookup_counts={
                "edit_graph_features": len(edit_graph_features),
            },
            index_id=resolved_index_id,
        )
        connection.execute(
            "INSERT INTO metadata(singleton, payload) VALUES (1, ?)",
            (_canonical_json_bytes(metadata).decode("utf-8"),),
        )
        connection.commit()
        connection.close()
        temporary.replace(destination)
        completed = True
    finally:
        try:
            connection.close()
        finally:
            if not completed and temporary.is_file():
                temporary.unlink()
    return {
        "schema_version": metadata["schema_version"],
        "storage_schema_version": metadata["storage_schema_version"],
        "precedent_scope": metadata["precedent_scope"],
        "index_id": metadata["index_id"],
        "row_count": metadata["row_count"],
        "lookup_counts": metadata["lookup_counts"],
        "path": str(destination),
    }


def load_sqlite_generic_index(path: str | Path) -> GenericReactionIndex:
    """Validate metadata and return a lazily materialized SQLite index."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    try:
        with _read_connection(source) as connection:
            result = connection.execute(
                "SELECT payload FROM metadata WHERE singleton = 1"
            ).fetchone()
    except sqlite3.DatabaseError as exc:
        raise ValueError("Not a valid SQLite recommendation index") from exc
    if result is None:
        raise ValueError("SQLite recommendation index metadata is missing")
    metadata = json.loads(result[0])
    if metadata.get("storage_schema_version") != SQLITE_INDEX_STORAGE_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported SQLite index storage schema; rebuild recommendation artifacts"
        )
    if metadata.get("row_encoding") != _ROW_ENCODING:
        raise ValueError("Unsupported SQLite precedent encoding; rebuild the index")
    if metadata.get("positions_encoding") != _POSITIONS_ENCODING:
        raise ValueError("Unsupported SQLite lookup encoding; rebuild the index")
    precedent_scope = _validate_index_metadata(metadata)
    row_count = int(metadata.get("row_count", -1))
    if row_count < 0:
        raise ValueError("Invalid SQLite recommendation index row count")
    lookup_counts = metadata.get("lookup_counts") or {}
    auxiliary_lookup_counts = metadata.get("auxiliary_lookup_counts") or {}
    map_names = (
        "exact",
        "handles",
        "transformations",
        "bond_edits",
        "environments",
        "facet_exact",
        "facet_attachments",
        "core_exact",
        "core_typed",
        "core_shapes",
        "core_centers",
        "environment_features",
        "departing_fragments",
        "fallback_features",
        "partial_transformations",
        "families",
    )
    if set(lookup_counts) != set(map_names):
        raise ValueError("SQLite recommendation lookup metadata is incomplete")
    if set(auxiliary_lookup_counts) != {"edit_graph_features"}:
        raise ValueError("SQLite auxiliary lookup metadata is incomplete")
    maps = {
        name: SQLiteLookupMapping(source, name, int(lookup_counts[name]))
        for name in map_names
    }
    return GenericReactionIndex(
        rows=SQLiteReactionRows(source, row_count),
        exact=maps["exact"],
        handles=maps["handles"],
        transformations=maps["transformations"],
        bond_edits=maps["bond_edits"],
        environments=maps["environments"],
        facet_exact=maps["facet_exact"],
        facet_attachments=maps["facet_attachments"],
        core_exact=maps["core_exact"],
        core_typed=maps["core_typed"],
        core_shapes=maps["core_shapes"],
        core_centers=maps["core_centers"],
        environment_features=maps["environment_features"],
        departing_fragments=maps["departing_fragments"],
        fallback_features=maps["fallback_features"],
        partial_transformations=maps["partial_transformations"],
        families=maps["families"],
        reaction_signature_schema_version=str(
            metadata["reaction_signature_schema_version"]
        ),
        reaction_core_schema_version=str(metadata["reaction_core_schema_version"]),
        reaction_core_algorithm_version=str(
            metadata["reaction_core_algorithm_version"]
        ),
        taxonomy_definition_versions=tuple(
            sorted(
                (str(key), str(value))
                for key, value in metadata["taxonomy_definition_versions"].items()
            )
        ),
        fallback_descriptor_schema_version=str(
            metadata["fallback_descriptor_schema_version"]
        ),
        fallback_definition_versions=tuple(
            sorted(
                (str(key), str(value))
                for key, value in metadata["fallback_definition_versions"].items()
            )
        ),
        record_schema_versions=tuple(metadata["record_schema_versions"]),
        converter_definition_versions=tuple(
            metadata["converter_definition_versions"]
        ),
        precedent_scope=precedent_scope,
        core_eligibility_definition_version=str(
            metadata["core_eligibility_definition_version"]
        ),
    )


def validate_sqlite_generic_index(path: str | Path) -> Dict[str, Any]:
    """Check SQLite structure, counts, compressed values, and row bounds."""
    source = Path(path)
    index = load_sqlite_generic_index(source)
    issues = []
    metadata_counts = {
        name: len(mapping) for name, mapping in _index_maps(index).items()
    }
    with _read_connection(source) as connection:
        metadata = json.loads(
            connection.execute(
                "SELECT payload FROM metadata WHERE singleton = 1"
            ).fetchone()[0]
        )
        auxiliary_counts = {
            str(name): int(count)
            for name, count in (
                metadata.get("auxiliary_lookup_counts") or {}
            ).items()
        }
        quick_check = str(connection.execute("PRAGMA quick_check").fetchone()[0])
        actual_row_count = int(
            connection.execute("SELECT COUNT(*) FROM precedents").fetchone()[0]
        )
        actual_lookup_counts = {
            str(kind): int(count)
            for kind, count in connection.execute(
                "SELECT lookup_kind, COUNT(*) FROM lookups GROUP BY lookup_kind"
            )
        }
        if quick_check != "ok":
            issues.append(f"sqlite_quick_check:{quick_check}")
        if actual_row_count != len(index.rows):
            issues.append("precedent_row_count_mismatch")
        for name, expected in metadata_counts.items():
            if actual_lookup_counts.get(name, 0) != expected:
                issues.append(f"lookup_count_mismatch:{name}")
        for name, expected in auxiliary_counts.items():
            if actual_lookup_counts.get(name, 0) != expected:
                issues.append(f"auxiliary_lookup_count_mismatch:{name}")
        for kind, key, payload in connection.execute(
            "SELECT lookup_kind, lookup_key, positions FROM lookups"
        ):
            try:
                positions = tuple(int(value) for value in _decompress_json(payload))
            except (TypeError, ValueError, zlib.error, json.JSONDecodeError):
                issues.append(f"lookup_payload_invalid:{kind}:{key}")
                continue
            if positions != tuple(sorted(set(positions))):
                issues.append(f"lookup_positions_not_sorted_unique:{kind}:{key}")
            if positions and (positions[0] < 0 or positions[-1] >= len(index.rows)):
                issues.append(f"lookup_position_out_of_bounds:{kind}:{key}")
    row_decode_count = 0
    try:
        for row_decode_count, _ in enumerate(index.rows, start=1):
            pass
    except (KeyError, TypeError, ValueError, zlib.error, json.JSONDecodeError):
        issues.append("precedent_payload_invalid")
    if row_decode_count != len(index.rows):
        issues.append("precedent_decode_count_mismatch")
    return {
        "schema_version": "1.0",
        "artifact_type": "sqlite_generic_index_integrity",
        "path": str(source),
        "valid": not issues,
        "issues": sorted(set(issues)),
        "index_schema_version": GENERIC_INDEX_SCHEMA_VERSION,
        "row_count": len(index.rows),
        "file_size_bytes": source.stat().st_size,
        "key_counts": metadata_counts,
    }


__all__ = [
    "SQLITE_INDEX_STORAGE_SCHEMA_VERSION",
    "SQLiteIndexBuildCancelled",
    "SQLiteLookupMapping",
    "SQLiteReactionRows",
    "build_sqlite_generic_index",
    "load_sqlite_generic_index",
    "save_sqlite_generic_index",
    "validate_sqlite_generic_index",
]
