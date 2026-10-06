"""Offline, immutable product-fragment index over canonical observations.

One SQLite artifact contains the RDKit library and metadata, so workspace content
pinning covers the entire index without a mutable companion-file dependency.
"""

from __future__ import annotations

from collections import Counter
from contextlib import closing
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
from time import monotonic
from typing import Any, Callable
import zlib
import zstandard

from rdkit import Chem, rdBase
from rdkit.Chem import rdSubstructLibrary

from reactive_taxonomy.fragment_search import (
    fragment_search_policy, indexed_product, project_fragment_evidence,
)
from .corpus_io import canonical_source_files, file_sha256, iter_canonical_records

import sqlite3

SCHEMA_VERSION = "fragment_precedent_index.v3"
LIBRARY_CHUNK_BYTES = 64 * 1024 * 1024
_RECORD_FIELDS = (
    "observation_id", "reaction_id", "reaction_smiles", "canonical_reaction_smiles",
    "reference_id", "reference_identity", "source_dataset", "source_path", "source_row_number",
    "admission_tier", "admission_reasons", "chemistry_status", "evidence_quality",
    "yield_pct", "temperature_c", "time_h", "conditions", "resolved_recipe", "synthesis_protocol",
    "source",
)


def pack(value: Any) -> bytes:
    """Compress JSON evidence without executable serialization."""
    return zlib.compress(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def unpack(value: bytes) -> Any:
    """Read JSON evidence stored by the index builder."""
    return json.loads(zlib.decompress(value))


def _write_serialized_library(connection: sqlite3.Connection, payload: bytes,
                              *, chunk_bytes: int = LIBRARY_CHUNK_BYTES) -> dict[str, Any]:
    """Store bounded compressed chunks, including incompressible large libraries."""
    encoder = zstandard.ZstdCompressor(level=6)
    connection.execute("DELETE FROM library")
    count = 0
    for ordinal, start in enumerate(range(0, len(payload), chunk_bytes)):
        chunk = payload[start:start + chunk_bytes]
        connection.execute("INSERT INTO library VALUES (?,?,?,?)", (
            ordinal, encoder.compress(chunk), len(chunk), hashlib.sha256(chunk).hexdigest()))
        count += 1
    return {"encoding": "zstandard_chunks.v1", "chunk_count": count,
            "serialized_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


def _read_serialized_library(connection: sqlite3.Connection, storage: dict[str, Any]) -> bytes:
    """Reject missing, reordered, corrupt, or incorrectly sized library chunks."""
    if storage.get("encoding") != "zstandard_chunks.v1":
        raise ValueError("Unsupported fragment library encoding")
    chunks = []
    decoder = zstandard.ZstdDecompressor()
    for expected, (ordinal, compressed, size, digest) in enumerate(connection.execute(
        "SELECT ordinal,payload,raw_size,sha256 FROM library ORDER BY ordinal"
    )):
        if ordinal != expected or not 0 < size <= LIBRARY_CHUNK_BYTES:
            raise ValueError("Fragment library chunk order or size is invalid")
        try:
            chunk = decoder.decompress(compressed, max_output_size=size)
        except zstandard.ZstdError as exc:
            raise ValueError("Fragment library chunk is corrupt") from exc
        if len(chunk) != size or hashlib.sha256(chunk).hexdigest() != digest:
            raise ValueError("Fragment library chunk checksum mismatch")
        chunks.append(chunk)
    payload = b"".join(chunks)
    if (len(chunks) != storage.get("chunk_count") or len(payload) != storage.get("serialized_bytes")
            or hashlib.sha256(payload).hexdigest() != storage.get("sha256")):
        raise ValueError("Fragment library is incomplete or its checksum mismatches")
    return payload


def load_fragment_library(connection: sqlite3.Connection, manifest: dict[str, Any]) -> Any:
    """Load the integrity-checked, chunked RDKit library from a completed index."""
    return rdSubstructLibrary.SubstructLibrary(
        _read_serialized_library(connection, manifest["library_storage"]))


def _finish_library(connection: sqlite3.Connection, library: Any,
                    manifest: dict[str, Any]) -> dict[str, Any]:
    """Finish packing without losing the completed source scan on a packing error."""
    completed = {k: v for k, v in manifest.items() if k not in {"scan_complete", "requested_max_records"}}
    completed["library_storage"] = _write_serialized_library(connection, library.Serialize())
    completed["build_complete"] = True
    completed["index_id"] = "FPI1:" + hashlib.sha256(json.dumps(completed, sort_keys=True).encode()).hexdigest()
    connection.execute("UPDATE metadata SET payload=?", (json.dumps(completed, sort_keys=True),))
    connection.commit()
    return completed


def build_fragment_index(
    source: str | Path, destination: str | Path, *, procedure_catalog: str | Path | None = None,
    evidence_catalog: str | Path | None = None,
    max_records: int | None = None, progress: Callable[[dict[str, Any]], None] | None = None,
    resume: bool = True,
) -> dict[str, Any]:
    """Build atomically; optional prefix pilots are explicitly marked in the manifest."""
    if max_records is not None and (type(max_records) is not int or max_records < 1):
        raise ValueError("max_records must be positive or omitted for full coverage")
    source, destination = Path(source).resolve(), Path(destination).resolve()
    sources = canonical_source_files(source, strict=True)
    inputs = list(dict.fromkeys((source, *sources, *([Path(procedure_catalog).resolve()] if procedure_catalog else []),
                                *([Path(evidence_catalog).resolve()] if evidence_catalog else []))))
    if destination in inputs:
        raise ValueError("Index destination cannot overwrite a source")
    identities = []
    for path in inputs:
        before = path.stat()
        digest = file_sha256(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("Canonical source changed during fingerprinting")
        identities.append({"path": str(path), "sha256": digest, "size_bytes": after.st_size,
                           "mtime_ns": after.st_mtime_ns})
    policy = fragment_search_policy()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".building")
    if temporary.exists() and resume:
        with closing(sqlite3.connect(temporary)) as connection:
            metadata_exists = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='metadata'").fetchone()
            pending = connection.execute("SELECT payload FROM metadata").fetchone() if metadata_exists else None
            pending = json.loads(pending[0]) if pending else None
            if pending and pending.get("scan_complete"):
                if (pending.get("sources") != identities or pending.get("policy") != policy
                    or pending.get("schema_version") != SCHEMA_VERSION
                    or pending.get("rdkit_version") != rdBase.rdkitVersion
                    or pending.get("requested_max_records") != max_records):
                    raise ValueError("Fragment staging inputs changed; rebuild with resume=False")
                library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.CachedMolHolder(),
                                                             rdSubstructLibrary.PatternHolder())
                if progress:
                    progress({"stage": "resume_library_packing", **pending["counts"]})
                for identity, smiles in connection.execute("SELECT id,smiles FROM products ORDER BY id"):
                    molecule = Chem.MolFromSmiles(smiles)
                    if molecule is None or library.AddMol(molecule) != identity:
                        raise ValueError("Fragment staging molecule library is inconsistent")
                report = _finish_library(connection, library, pending)
            else:
                report = None
        if report is not None:
            os.replace(temporary, destination)
            return report
    # A partial scan has no complete checkpoint; only this builder's owned
    # staging file is replaced. Published outputs remain untouched.
    temporary.unlink(missing_ok=True)
    library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.CachedMolHolder(),
                                                 rdSubstructLibrary.PatternHolder())
    counts: Counter = Counter()
    products: dict[str, int] = {}
    started = monotonic()

    # Repeated experiments commonly share identical structures and correspondence.
    # Cache only inputs actually consumed by the projection; never merge observations.
    @lru_cache(maxsize=2048)
    def evidence_for(raw: str, observation_json: str) -> list[dict[str, Any]]:
        return project_fragment_evidence(raw, json.loads(observation_json))

    normalize_product = lru_cache(maxsize=4096)(indexed_product)
    with closing(sqlite3.connect(temporary)) as connection:
        connection.executescript("""
            CREATE TABLE metadata (payload TEXT NOT NULL);
            CREATE TABLE library (ordinal INTEGER PRIMARY KEY, payload BLOB NOT NULL,
                raw_size INTEGER NOT NULL, sha256 TEXT NOT NULL);
            CREATE TABLE products (id INTEGER PRIMARY KEY, smiles TEXT UNIQUE NOT NULL);
            CREATE TABLE observations (id TEXT PRIMARY KEY, reaction_id TEXT, reference_id TEXT, payload BLOB);
            CREATE TABLE links (product_id INTEGER, observation_id TEXT, component_index INTEGER,
                atom_order TEXT, evidence BLOB, side TEXT NOT NULL,
                PRIMARY KEY(observation_id, side, component_index));
            CREATE INDEX links_product ON links(product_id);
            CREATE INDEX links_side_product ON links(side,product_id);
            CREATE TABLE procedures (observation_id TEXT, reaction_id TEXT, payload BLOB);
            CREATE INDEX procedures_observation ON procedures(observation_id);
            CREATE INDEX procedures_reaction ON procedures(reaction_id);
        """)
        prefix_limited = False
        for record in iter_canonical_records(source, strict=True):
            if max_records is not None and counts["source_observations"] >= max_records:
                prefix_limited = True
                break
            counts["source_observations"] += 1
            identity = record.get("observation_id")
            if not isinstance(identity, str) or not identity:
                raise ValueError("Canonical observation_id is required; reconvert this source")
            raw = str(record.get("reaction_smiles") or "")
            parts = raw.split(">")
            if len(parts) != 3:
                counts["invalid_reactions"] += 1
                continue
            observation = record.get("reaction_observation") or {}
            projection_inputs = {key: observation.get(key) for key in (
                "valid", "input_reaction_smiles", "evidence_quality", "warnings",
            )}
            projection_inputs["edit_hypotheses"] = bool(observation.get("edit_hypotheses") or record.get("reaction_edit_hypotheses"))
            projection_inputs["evidence_candidates"] = [
                {"status": c.get("status")} for c in
                (observation.get("evidence_candidates") or record.get("reaction_evidence_candidates") or ())]
            projected = evidence_for(raw, json.dumps(projection_inputs, sort_keys=True))
            fields = ("observation_id", "reaction_id", "reference_id", "reaction_smiles",
                      "admission_tier", "admission_reasons", "evidence_quality") if evidence_catalog else _RECORD_FIELDS
            details = {key: record.get(key) for key in fields}
            details["discovery_projection_sha256"] = hashlib.sha256(json.dumps(details, sort_keys=True).encode()).hexdigest()
            details["warnings"] = list((record.get("reaction_observation") or {}).get("warnings") or ())
            connection.execute("INSERT INTO observations VALUES (?,?,?,?)", (
                identity, record.get("reaction_id", ""), record.get("reference_id", ""), pack(details)))
            for side, ci, smiles in (
                (side, ci, smiles)
                for side, index in (("product", 2), ("reactant", 0))
                for ci, smiles in enumerate(parts[index].split(".")) if smiles
            ):
                try:
                    canonical, order = normalize_product(smiles)
                except ValueError:
                    counts[f"invalid_{side}_components"] += 1
                    continue
                product_id = products.get(canonical)
                if product_id is None:
                    product_id = library.AddMol(Chem.MolFromSmiles(canonical))
                    products[canonical] = product_id
                    connection.execute("INSERT INTO products VALUES (?,?)", (product_id, canonical))
                evidence = (projected[ci] if ci < len(projected) else {"evidence_status": "unresolved"}) if side == "product" else {
                    "evidence_status": "reported_reactant_occurrence",
                    "warnings": list(observation.get("warnings") or ()),
                }
                connection.execute("INSERT INTO links VALUES (?,?,?,?,?,?)", (
                    product_id, identity, ci, json.dumps(order), pack(evidence), side))
                counts[f"indexed_{side}_components"] += 1
                counts["evidence_" + evidence["evidence_status"]] += 1
            if counts["source_observations"] % 1000 == 0:
                connection.commit()
                if progress:
                    progress({"stage": "index_observations", **counts,
                              "elapsed_seconds": round(monotonic() - started, 3)})
        if procedure_catalog:
            for record in iter_canonical_records(procedure_catalog):
                connection.execute("INSERT INTO procedures VALUES (?,?,?)", (
                    record.get("observation_id") or "", record.get("reaction_id") or "", pack(record)))
                counts["procedures"] += 1
        for item in identities:
            stat = Path(item["path"]).stat()
            if (stat.st_size, stat.st_mtime_ns) != (item["size_bytes"], item["mtime_ns"]):
                raise ValueError("Canonical source changed during index build")
        counts["distinct_products"] = len(products)
        manifest = {"schema_version": SCHEMA_VERSION, "definition_version": policy["definition_version"],
                    "policy": policy, "rdkit_version": rdBase.rdkitVersion,
                    "build_complete": False, "source_coverage_complete": not prefix_limited,
                    "source_scope": "prefix_pilot" if prefix_limited else "full_source",
                    "counts": dict(counts), "sources": identities,
                    "record_scope": "discovery_fields_not_complete_canonical_record",
                    "correspondence_scope": "validated_supplied_maps_only"}
        manifest["search_sides"] = ["product", "reactant"]
        manifest["counts"]["distinct_molecules"] = len(products)
        manifest["counts"]["distinct_products"] = connection.execute(
            "SELECT count(DISTINCT product_id) FROM links WHERE side='product'").fetchone()[0]
        if evidence_catalog:
            catalog_path = Path(evidence_catalog).resolve()
            with closing(sqlite3.connect(catalog_path.as_uri() + "?mode=ro", uri=True)) as catalog_db:
                catalog_metadata = json.loads(catalog_db.execute("SELECT payload FROM metadata").fetchone()[0])
            manifest["evidence_catalog"] = {
                "relative_path": os.path.relpath(catalog_path, destination.parent),
                "sha256": next(i["sha256"] for i in identities if i["path"] == str(catalog_path)),
                "size_bytes": catalog_path.stat().st_size,
                "catalog_id": catalog_metadata["catalog_id"],
            }
            manifest["counts"]["procedures"] = catalog_metadata["counts"]["procedures"]
        manifest.update(scan_complete=True, requested_max_records=max_records)
        connection.execute("INSERT INTO metadata VALUES (?)", (json.dumps(manifest, sort_keys=True),))
        connection.commit()
        if progress:
            progress({"stage": "library_packing", **counts})
        manifest = _finish_library(connection, library, manifest)
    os.replace(temporary, destination)
    return manifest



def open_fragment_index(path: str | Path) -> tuple[sqlite3.Connection, dict[str, Any]]:
    """Open a completed compatible immutable artifact; never rebuild at search time."""
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError("Fragment index unavailable; build with python -m condition_recommender.fragment_search build")
    connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
    try:
        manifest = json.loads(connection.execute("SELECT payload FROM metadata").fetchone()[0])
        if (manifest.get("schema_version") != SCHEMA_VERSION or not manifest.get("build_complete")
                or manifest.get("rdkit_version") != rdBase.rdkitVersion
                or manifest.get("policy") != fragment_search_policy()):
            raise ValueError("Fragment index is incompatible or incomplete; rebuild offline")
        dependency = manifest.get("evidence_catalog")
        if dependency:
            catalog_path = (path.parent / dependency["relative_path"]).resolve()
            if not catalog_path.is_file() or catalog_path.stat().st_size != dependency["size_bytes"]:
                raise ValueError("Fragment evidence catalog is missing or incompatible")
            with closing(sqlite3.connect(catalog_path.as_uri() + "?mode=ro", uri=True)) as catalog_db:
                catalog_metadata = json.loads(catalog_db.execute("SELECT payload FROM metadata").fetchone()[0])
            if not catalog_metadata.get("build_complete") or catalog_metadata.get("catalog_id") != dependency["catalog_id"]:
                raise ValueError("Fragment evidence catalog binding mismatch")
        return connection, manifest
    except BaseException:
        connection.close()
        raise
