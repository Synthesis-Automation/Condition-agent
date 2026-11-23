"""Indexed access to shared canonical evidence without scanning gzip catalogs."""

from __future__ import annotations

from collections import Counter, OrderedDict
from collections.abc import Mapping, Iterator
from contextlib import closing
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sqlite3
from typing import Any, Iterable, Callable
import zstandard

from reactive_taxonomy.shared_reaction_core import build_shared_reaction_core, shared_core_definition_hash

from .conversion.atomic import atomic_output_path
from .conversion.identities import canonical_reaction_identity
from .corpus_io import canonical_source_files, file_sha256
from .record_storage import iter_record_events, json_bytes

CATALOG_SCHEMA_VERSION = "processed_evidence_catalog.v3"
_LEGACY_SCHEMA = "processed_evidence_catalog.v2"
EXPERIMENT_FIELDS = (
    "observation_id", "reaction_id", "reaction_smiles", "canonical_reaction_id",
    "reference_id", "reference_identity", "source_dataset", "source_path", "source_row_number",
    "source", "resolved_recipe", "resolved_recipe_id", "outcome_status", "yield_pct",
    "temperature_c", "time_h", "conditions",
    "index_eligibility", "admission_tier", "admission_reasons", "chemistry_status",
)


@dataclass(frozen=True)
class ExperimentLookup:
    """Source records in a bounded exact or graph-local experimental search."""

    records: tuple[dict[str, Any], ...]
    truncated: bool
    candidates_examined: int
    lookup_complete: bool
    scope: str


def build_processed_catalog(source: str | Path, destination: str | Path, *,
                            include_local_cores: bool = True,
                            progress: Callable[[dict[str, Any]], None] | None = None) -> dict[str, Any]:
    """Materialize canonical shared objects and indexed observation metadata once."""
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    shards = canonical_source_files(source, strict=True)
    if output.resolve() == Path(source).resolve() or output.resolve() in shards:
        raise ValueError("Evidence catalog destination cannot overwrite a source")
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

    # Repeated wells share a graph but each has a distinct source observation.
    @lru_cache(maxsize=1024)
    def exact_key(smiles: str) -> str | None:
        identity = canonical_reaction_identity(smiles)
        return identity.reaction_id if identity is not None else None

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
                CREATE TABLE experiment_lookup(kind TEXT NOT NULL, key TEXT NOT NULL,
                    observation_id TEXT NOT NULL,
                    PRIMARY KEY(kind, key, observation_id)) WITHOUT ROWID;
                CREATE TABLE reference_records(id TEXT PRIMARY KEY,payload BLOB NOT NULL) WITHOUT ROWID;
                CREATE TABLE recipe_records(id TEXT PRIMARY KEY,payload BLOB NOT NULL) WITHOUT ROWID;
            """)
            db.execute("INSERT INTO compression_dictionary VALUES (?)", (dictionary,))
            for shard in shards:
                local: dict[str, Any] = {}
                projections: OrderedDict[tuple, tuple[tuple[str, str], ...]] = OrderedDict()

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

                def local_keys(row: Mapping[str, Any], smiles: str) -> tuple[tuple[str, str], ...]:
                    # Reuse the verified conversion evidence. Repeated wells normally
                    # share object IDs, so project their graph only once per shard.
                    if not row.get("reaction_signature") or not row.get("reaction_core"):
                        return ()
                    if restore(row.get("reaction_completeness") or {}).get("status") != "verified":
                        return ()
                    inputs = (row["reaction_signature"], row["reaction_core"])
                    key = (smiles, *(item.get("$reaction_object") or hashlib.sha256(json_bytes(item)).digest()
                                     for item in inputs))
                    if key in projections:
                        projections.move_to_end(key)
                        return projections[key]
                    projection = build_shared_reaction_core(smiles, *(restore(item) for item in inputs))
                    levels = tuple((level.level, level.key) for level in projection.levels)
                    projections[key] = levels
                    if len(projections) > 1024:
                        projections.popitem(last=False)
                    return levels

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
                    smiles = str(row.get("reaction_smiles") or "")
                    identity = exact_key(smiles)
                    keys = [("whole_reaction", identity)] if identity else []
                    if identity and include_local_cores:
                        levels = local_keys(row, smiles)
                        keys.extend(levels)
                        counts["experiments_with_local_core"] += bool(levels)
                    db.executemany("INSERT INTO experiment_lookup VALUES (?,?,?)",
                                   ((kind, key, oid) for kind, key in keys))
                    counts["experiments_with_graph"] += bool(keys)
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
            local_hash = shared_core_definition_hash() if include_local_cores else None
            catalog_id = hashlib.sha256(json_bytes(
                [file_sha256(Path(source)), CATALOG_SCHEMA_VERSION, local_hash]
            )).hexdigest()
            report = {"schema_version": CATALOG_SCHEMA_VERSION, "build_complete": True,
                      "catalog_id": catalog_id,
                      "experiment_local_core_definition": local_hash,
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
        self._file_identity = self._identity()

    def _identity(self) -> tuple[int, int, int]:
        stat = self.path.stat()
        return stat.st_ino, stat.st_size, stat.st_mtime_ns

    def _connect(self) -> sqlite3.Connection:
        db = sqlite3.connect(self.path.as_uri() + "?mode=ro&immutable=1", uri=True)
        try:
            if self._identity() != self._file_identity:
                raise ValueError("Evidence catalog changed; reopen it before reading")
            return db
        except BaseException:
            db.close()
            raise

    @lru_cache(maxsize=1)
    def _metadata(self) -> dict[str, Any]:
        with closing(self._connect()) as db:
            try:
                row = db.execute("SELECT payload FROM metadata").fetchone()
                metadata = json.loads(row[0]) if row else None
            except (sqlite3.DatabaseError, ValueError) as exc:
                raise ValueError("File is not a readable evidence catalog") from exc
        if (not isinstance(metadata, dict)
                or metadata.get("schema_version") not in (CATALOG_SCHEMA_VERSION, _LEGACY_SCHEMA)
                or not metadata.get("build_complete")):
            raise ValueError("Evidence catalog is incompatible or incomplete")
        return metadata

    @property
    def supports_experiment_analogues(self) -> bool:
        """Whether the artifact has a graph-local experimental lookup.

        A changed chemistry definition requires rebuilding the derived keys.
        """
        metadata = self._metadata()
        definition = metadata.get("experiment_local_core_definition")
        if metadata["schema_version"] == _LEGACY_SCHEMA or definition is None:
            return False
        if definition != shared_core_definition_hash():
            raise ValueError("Experiment local-core definition mismatch; rebuild the experiment catalog")
        return True

    @lru_cache(maxsize=1)
    def _dictionary(self) -> bytes:
        metadata = self._metadata()
        with closing(self._connect()) as db:
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

    def experiment_candidates(
        self, reaction_smiles: str, *, fragment_index: str | Path | None = None,
        limit: int = 128, include_analogues: bool = False,
    ) -> ExperimentLookup:
        """Find observations of the intended graph regardless of recipe admission.

        Latest catalogs also index verified local reaction edits. Candidate
        generation precedes structural qualification and recipe comparison. A
        yield measurement of zero never establishes formation of the target.
        """
        if type(limit) is not int or not 1 <= limit <= 512:
            raise ValueError("experiment candidate limit must be between 1 and 512")
        metadata = self._metadata()
        legacy = metadata["schema_version"] == _LEGACY_SCHEMA
        if include_analogues:
            if legacy and fragment_index is None:
                raise ValueError("Legacy experiment catalog requires its bound fragment index")
            if not legacy and not self.supports_experiment_analogues:
                raise ValueError("Local experiment lookup unavailable; rebuild the experiment catalog")
        identity = canonical_reaction_identity(reaction_smiles)
        if identity is None:
            return ExperimentLookup((), False, 0, False, "invalid_reaction_graph")
        if legacy:
            if fragment_index is None:
                return ExperimentLookup((), False, 0, False, "exact_graph_index_unavailable")
            return self._legacy_experiments(identity.reaction_id, identity.products, fragment_index, limit)

        keys = [("whole_reaction", identity.reaction_id)]
        scope = "exact_graph"
        if include_analogues:
            from .reaction_precedent_support import reaction_evidence_core

            core = reaction_evidence_core(reaction_smiles)
            if core is not None and core.levels:
                keys.extend((level.level, level.key) for level in core.levels)
                scope = "exact_and_indexed_local_cores"
        found: dict[str, None] = {}
        with closing(self._connect()) as db:
            for kind, key in keys:
                cursor = db.execute(
                    "SELECT observation_id FROM experiment_lookup WHERE kind=? AND key=? "
                    "ORDER BY observation_id LIMIT ?", (kind, key, limit + 1),
                )
                for (obs_id,) in cursor:
                    found[obs_id] = None
                    if len(found) > limit:
                        break
                if len(found) > limit:
                    break
        truncated = len(found) > limit
        records = tuple(self.observation(i, EXPERIMENT_FIELDS) for i in list(found)[:limit])
        return ExperimentLookup(records, truncated, len(records), not truncated, scope)

    def _legacy_experiments(
        self, reaction_id: str, products: tuple[str, ...], fragment_index: str | Path, limit: int,
    ) -> ExperimentLookup:
        """Use a catalog-bound product index and verify both molecular sides."""
        from reactive_taxonomy.fragment_search import indexed_product
        from .fragment_index import open_fragment_index, unpack

        fragment_db, fragment = open_fragment_index(fragment_index)
        with closing(fragment_db):
            binding = fragment.get("evidence_catalog") or {}
            if binding.get("catalog_id") != self._metadata()["catalog_id"]:
                raise ValueError("Fragment index is not bound to this experiment catalog")
            # The full product multiset, stereo and both precursor identity and
            # multiplicity are verified below, including for multi-product steps.
            product_key, _ = indexed_product(products[0])
            scan_limit = min(512, 4 * limit)
            occurrences = fragment_db.execute(
                "SELECT observation_id FROM links WHERE side='product' AND product_id="
                "(SELECT id FROM products WHERE smiles=?) ORDER BY rowid LIMIT ?",
                (product_key, scan_limit + 1),
            ).fetchall()
            capped = len(occurrences) > scan_limit
            visited = set()
            matched = []
            for (obs_id,) in occurrences[:scan_limit]:
                if obs_id in visited:
                    continue
                visited.add(obs_id)
                packed = fragment_db.execute("SELECT payload FROM observations WHERE id=?", (obs_id,)).fetchone()
                source = canonical_reaction_identity(unpack(packed[0])["reaction_smiles"])
                if source is not None and source.reaction_id == reaction_id:
                    matched.append(obs_id)
                    if len(matched) > limit:
                        capped = True
                        break
        records = tuple(self.observation(i, EXPERIMENT_FIELDS) for i in matched[:limit])
        return ExperimentLookup(records, capped, len(visited),
                                not capped and bool(fragment.get("source_coverage_complete")),
                                "bound_fragment_exact_graph")

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
