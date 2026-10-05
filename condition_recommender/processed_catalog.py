"""Indexed access to shared canonical evidence without scanning gzip catalogs."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Iterator
from contextlib import closing
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sqlite3
from typing import Any, Iterable, Callable
import zstandard

from .conversion.atomic import atomic_output_path
from .corpus_io import canonical_source_files, file_sha256
from .record_storage import iter_record_events, json_bytes

CATALOG_SCHEMA_VERSION = "processed_evidence_catalog.v2"


def build_processed_catalog(source: str | Path, destination: str | Path, *,
                            progress: Callable[[dict[str, Any]], None] | None = None) -> dict[str, Any]:
    """Materialize canonical shared objects and indexed observation metadata once."""
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    shards = canonical_source_files(source, strict=True)
    samples = []
    for shard in shards:
        for event in iter_record_events(shard):
            if "object_id" in event:
                samples.append(json_bytes(event["value"]))
            elif "record" in event:
                samples.append(json_bytes(event["record"]))
            if len(samples) >= 2048:
                break
        if len(samples) >= 2048:
            break
    # Compression changes representation only. Store its exact dictionary in
    # the artifact so all evidence remains self-contained and lossless.
    dictionary = (zstandard.train_dictionary(32768, samples).as_bytes()
                  if len(samples) >= 256 and sum(map(len, samples)) >= 65536 else b"")
    compressor = zstandard.ZstdCompressor(level=6,
                    dict_data=zstandard.ZstdCompressionDict(dictionary) if dictionary else None)

    def encode(value: Any) -> bytes:
        return compressor.compress(json_bytes(value))

    with atomic_output_path(output) as temporary:
        with closing(sqlite3.connect(temporary)) as db:
            db.executescript("""
                PRAGMA journal_mode=OFF;
                PRAGMA synchronous=OFF;
                PRAGMA temp_store=FILE;
                CREATE TABLE metadata(payload TEXT NOT NULL);
                CREATE TABLE compression_dictionary(payload BLOB NOT NULL);
                CREATE TABLE objects(id TEXT PRIMARY KEY,payload BLOB NOT NULL) WITHOUT ROWID;
                CREATE TABLE observations(id TEXT PRIMARY KEY,reaction_id TEXT NOT NULL,
                    reference_id TEXT NOT NULL,source_dataset TEXT NOT NULL,
                    has_procedure INTEGER NOT NULL,payload BLOB NOT NULL) WITHOUT ROWID;
                CREATE INDEX observation_reaction ON observations(reaction_id);
                CREATE INDEX observation_reference ON observations(reference_id);
                CREATE TABLE reference_records(id TEXT PRIMARY KEY,payload BLOB NOT NULL) WITHOUT ROWID;
                CREATE TABLE recipe_records(id TEXT PRIMARY KEY,payload BLOB NOT NULL) WITHOUT ROWID;
            """)
            db.execute("INSERT INTO compression_dictionary VALUES (?)", (dictionary,))
            for shard in shards:
                local: dict[str, Any] = {}

                def restore(value: Any) -> Any:
                    if isinstance(value, dict):
                        if set(value) == {"$reaction_object"}:
                            identity = value["$reaction_object"]
                            if identity not in local:
                                raise ValueError(f"Missing shard object {identity}")
                            return restore(local[identity])
                        return {k: restore(v) for k, v in value.items()}
                    if isinstance(value, list):
                        return [restore(v) for v in value]
                    return value

                for event in iter_record_events(shard):
                    if "reset_objects" in event:
                        local.clear()
                        continue
                    if "object_id" in event:
                        local[event["object_id"]] = event["value"]
                        db.execute("INSERT OR IGNORE INTO objects VALUES (?,?)",
                                   (event["object_id"], encode(event["value"])))
                        continue
                    row = event["record"]
                    oid = str(row.get("observation_id") or "")
                    if not oid:
                        raise ValueError("Canonical observation ID is required")
                    source_row = restore(row.get("source") or {})
                    db.execute("INSERT INTO observations VALUES (?,?,?,?,?,?)", (
                        oid, str(row.get("reaction_id") or ""), str(row.get("reference_id") or ""),
                        str(row.get("source_dataset") or ""),
                        int(bool(str(source_row.get("experimental_procedure") or "").strip())),
                        encode(row)))
                    counts["observations"] += 1
                    counts["procedures"] += int(bool(str(source_row.get("experimental_procedure") or "").strip()))
                    for table, id_field, value_field in (
                        ("reference_records", "reference_id", "reference_identity"),
                        ("recipe_records", "resolved_recipe_id", "resolved_recipe"),
                    ):
                        identity = row.get(id_field)
                        if identity and row.get(value_field):
                            db.execute(f"INSERT OR IGNORE INTO {table} VALUES (?,?)",
                                       (identity, encode(row[value_field])))
                db.commit()
                if progress:
                    progress({"observations": counts["observations"], "shard": shard.name})
            for key, table in (("shared_objects", "objects"), ("references", "reference_records"),
                               ("recipes", "recipe_records")):
                counts[key] = db.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            catalog_id = hashlib.sha256((file_sha256(Path(source)) + CATALOG_SCHEMA_VERSION).encode()).hexdigest()
            report = {"schema_version": CATALOG_SCHEMA_VERSION, "build_complete": True,
                      "catalog_id": catalog_id,
                      "payload_encoding": "zstandard_dictionary.v1",
                      "dictionary_sha256": hashlib.sha256(dictionary).hexdigest(),
                      "counts": dict(counts)}
            db.execute("INSERT INTO metadata VALUES (?)", (json.dumps(report, sort_keys=True),))
            db.commit()
    return report


class ProcessedCatalog:
    """Read immutable observation and evidence objects through indexed IDs."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path).resolve()

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path.as_uri() + "?mode=ro&immutable=1", uri=True)

    @lru_cache(maxsize=1)
    def _dictionary(self) -> bytes:
        with closing(self._connect()) as db:
            metadata = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
            if metadata.get("schema_version") != CATALOG_SCHEMA_VERSION or not metadata.get("build_complete"):
                raise ValueError("Evidence catalog is incompatible or incomplete")
            dictionary = db.execute("SELECT payload FROM compression_dictionary").fetchone()[0]
        if hashlib.sha256(dictionary).hexdigest() != metadata.get("dictionary_sha256"):
            raise ValueError("Evidence compression dictionary checksum mismatch")
        return dictionary

    def _decode(self, payload: bytes) -> Any:
        dictionary = self._dictionary()
        decoder = zstandard.ZstdDecompressor(dict_data=zstandard.ZstdCompressionDict(dictionary) if dictionary else None)
        return json.loads(decoder.decompress(payload))

    @lru_cache(maxsize=4096)
    def _object(self, identity: str) -> Any:
        with closing(self._connect()) as db:
            row = db.execute("SELECT payload FROM objects WHERE id=?", (identity,)).fetchone()
        if row is None:
            raise ValueError(f"Missing canonical evidence object: {identity}")
        return self._decode(row[0])

    def restore(self, value: Any) -> Any:
        """Hydrate fresh nested values without exposing shared mutable cache entries."""
        if isinstance(value, dict):
            if set(value) == {"$reaction_object"}:
                return self.restore(self._object(value["$reaction_object"]))
            return {k: self.restore(v) for k, v in value.items()}
        if isinstance(value, list):
            return [self.restore(v) for v in value]
        return value

    def observation(self, identity: str, fields: Iterable[str] | None = None) -> dict[str, Any]:
        """Fetch one observation, optionally hydrating only selected top-level fields."""
        with closing(self._connect()) as db:
            row = db.execute("SELECT payload FROM observations WHERE id=?", (identity,)).fetchone()
        if row is None:
            raise KeyError(identity)
        value = self._decode(row[0])
        if fields is not None:
            selected = set(fields) | {"observation_id", "reaction_id"}
            value = {k: v for k, v in value.items() if k in selected}
        return self.restore(value)

    def references(self, identities: Iterable[str]) -> dict[str, dict[str, Any]]:
        """Fetch exact source reference records without loading the entire catalog."""
        result = {}
        with closing(self._connect()) as db:
            for identity in dict.fromkeys(identities):
                row = db.execute("SELECT payload FROM reference_records WHERE id=?", (identity,)).fetchone()
                if row is not None:
                    result[identity] = self.restore(self._decode(row[0]))
        return result

    def observation_procedure(self, identity: str) -> dict[str, Any] | None:
        """Read the procedure from its exact observation; never join by a broad alias."""
        row = self.observation(identity, ("source", "reference_id", "source_dataset"))
        source = row.pop("source")
        text = str(source.get("experimental_procedure") or "").strip()
        if not text:
            return None
        return {**row, "procedure_text": text, "notes": source.get("notes"),
                "stages": source.get("stages"), "steps": source.get("steps")}

    def procedures(self, reaction_ids: Iterable[str], *, offset: int = 0,
                   limit: int = 100) -> dict[str, Any]:
        """Return a bounded page of source procedures with exact observation links."""
        identities = tuple(dict.fromkeys(reaction_ids))
        if not identities:
            return {"records": [], "total": 0, "offset": offset, "next_offset": None, "missing_reaction_ids": []}
        if len(identities) > 100 or offset < 0 or not 1 <= limit <= 100:
            raise ValueError("Procedure query exceeds supported bounds")
        marks = ",".join("?" for _ in identities)
        where = f"has_procedure=1 AND reaction_id IN ({marks})"
        with closing(self._connect()) as db:
            total = db.execute(f"SELECT count(*) FROM observations WHERE {where}", identities).fetchone()[0]
            rows = db.execute(f"SELECT id FROM observations WHERE {where} ORDER BY id LIMIT ? OFFSET ?",
                              (*identities, limit, offset)).fetchall()
            found = {row[0] for row in db.execute(f"SELECT DISTINCT reaction_id FROM observations WHERE {where}", identities)}
        result = []
        for (identity,) in rows:
            row = self.observation(identity, ("source", "source_dataset", "reference_id"))
            source = row.pop("source")
            result.append({**row, "procedure_text": source.get("experimental_procedure"),
                           "notes": source.get("notes"), "stages": source.get("stages"),
                           "steps": source.get("steps")})
        return {"records": result, "total": total, "offset": offset,
                "missing_reaction_ids": sorted(set(identities) - found),
                "next_offset": offset + limit if offset + limit < total else None}


class ReferenceCatalogView(Mapping[str, dict[str, Any]]):
    """Mapping interface backed by exact indexed reference lookups."""

    def __init__(self, path: str | Path) -> None:
        self.catalog = ProcessedCatalog(path)

    def __getitem__(self, identity: str) -> dict[str, Any]:
        result = self.catalog.references((identity,))
        if identity not in result:
            raise KeyError(identity)
        return result[identity]

    def __len__(self) -> int:
        with closing(self.catalog._connect()) as db:
            return db.execute("SELECT count(*) FROM reference_records").fetchone()[0]

    def __iter__(self) -> Iterator[str]:
        with closing(self.catalog._connect()) as db:
            for (identity,) in db.execute("SELECT id FROM reference_records ORDER BY id"):
                yield identity


class ProcedureCatalogView(Mapping[str, dict[str, Any]]):
    """Read exact observation or explicit reaction-level procedure keys lazily."""

    def __init__(self, path: str | Path) -> None:
        self.catalog = ProcessedCatalog(path)

    def __getitem__(self, key: str) -> dict[str, Any]:
        kind, identity = key.split(":", 1)
        if kind == "observation":
            result = self.catalog.observation_procedure(identity)
        elif kind == "reaction":
            result = next(iter(self.catalog.procedures((identity,), limit=1)["records"]), None)
        else:
            result = None
        if result is None:
            raise KeyError(key)
        return result

    def __len__(self) -> int:
        with closing(self.catalog._connect()) as db:
            return db.execute("SELECT count(*) FROM observations WHERE has_procedure=1").fetchone()[0]

    def __iter__(self) -> Iterator[str]:
        with closing(self.catalog._connect()) as db:
            for (identity,) in db.execute("SELECT id FROM observations WHERE has_procedure=1 ORDER BY id"):
                yield "observation:" + identity
