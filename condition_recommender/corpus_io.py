"""Shared streaming access to canonical observations and validated shard manifests."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator
from .record_storage import iter_record_shard


def file_sha256(path: Path) -> str:
    """Hash a source in bounded memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_source_files(path: str | Path, *, strict: bool = False) -> tuple[Path, ...]:
    """Resolve manifest leaves; strict builds require complete, checksum-valid shards."""
    def resolve(source: Path, ancestors: frozenset[Path]) -> list[Path]:
        source = source.resolve()
        if source in ancestors:
            raise ValueError(f"Cyclic canonical manifest: {source}")
        if source.suffix != ".json":
            if not source.is_file():
                raise FileNotFoundError(source)
            return [source]
        payload = json.loads(source.read_text("utf-8"))
        if strict and (payload.get("coverage_complete") is False or any(
            not item.get("coverage_complete") for item in payload.get("source_files", ())
        )):
            raise ValueError(f"Incomplete canonical source coverage: {source}")
        kind = payload.get("artifact_type")
        if payload.get("schema_version") in {"processed_reaction_active.v1", "processed_reaction_release.v1"}:
            from .processed_release import resolve_processed_release
            release = resolve_processed_release(source)
            return resolve(release.artifact("canonical_records"), ancestors | {source})
        if kind == "saved_recommendation_batch_manifest":
            entries = payload.get("source_manifests") or ()
            if not entries:
                raise ValueError(f"Saved batch has no converted-source references: {source}")
            paths = [source.parent / e["relative_path"] if e.get("relative_path") else Path(e["path"])
                     for e in entries]
        elif kind == "generic_sharded_conversion":
            entries = payload.get("shards") or ()
            if strict and (not entries or any(e.get("status") != "complete" for e in entries)):
                raise ValueError(f"Incomplete canonical manifest: {source}")
            paths = []
            for entry in entries:
                if entry.get("status") != "complete":
                    continue
                target = source.parent / entry["output_path"]
                if strict and (not target.is_file() or file_sha256(target) != entry.get("output_sha256")):
                    raise ValueError(f"Canonical shard checksum failed: {target}")
                paths.append(target)
        else:
            raise ValueError(f"Not a sharded conversion manifest: {source}")
        return [leaf for target in paths for leaf in resolve(target, ancestors | {source})]
    paths = resolve(Path(path), frozenset())
    if strict and len(set(paths)) != len(paths):
        raise ValueError("Canonical manifest repeats a shard")
    return tuple(paths)


def iter_canonical_records(path: str | Path, *, strict: bool = False) -> Iterator[dict[str, Any]]:
    """Stream canonical JSONL/gzip or manifests without converting source chemistry."""
    for source in canonical_source_files(path, strict=strict):
        yield from iter_record_shard(source)
