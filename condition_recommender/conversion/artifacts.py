"""Build scalable canonical and indexed recommendation artifacts."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, Mapping, Optional

from ..sqlite_indexing import (
    SQLiteIndexBuildCancelled,
    build_sqlite_generic_index,
)
from .concise_review import iter_canonical_records
from .sharded import (
    SHARD_MANIFEST_SCHEMA_VERSION,
    SHARDED_CONVERSION_DEFINITION_VERSION,
    ShardedConversionCancelled,
    ShardedConversionProgress,
    convert_datasets_sharded,
    iter_gzip_jsonl,
    validate_sharded_conversion,
    write_conversion_catalogs,
)
from .input_schema import ConversionDatasetInput, discover_conversion_datasets

RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION = "2.3"
SAVED_BATCH_WORKFLOW_SCHEMA_VERSION = "1.3"
SAVED_BATCHES_DIRNAME = "batches"
CONVERTED_SOURCES_DIRNAME = "converted_sources"
SAVED_BATCH_MANIFEST_ARTIFACT_TYPE = "saved_recommendation_batch_manifest"
COMBINED_RECORDS_FILENAME = "combined_records.jsonl.gz"
COMBINED_BATCH_MANIFEST_FILENAME = "combined_batch_manifest.json"
RECOMMENDATION_LIBRARY_MODES = ("full", "compact")


def recommendation_library_mode_dir(
    library_root: str | Path,
    mode: str,
) -> Path:
    """Return the isolated artifact directory for a library build mode.

    Existing root-level libraries predate the mode layout and remain the Full
    location until explicitly migrated. This avoids requiring an expensive Full
    rebuild merely to introduce a Compact development library.
    """
    normalized = mode.strip().casefold()
    if normalized not in RECOMMENDATION_LIBRARY_MODES:
        raise ValueError(f"Unsupported recommendation library mode: {mode}")
    root = Path(library_root)
    mode_dir = root / normalized
    legacy_full_markers = (
        root / "shard_manifest.json",
        root / SAVED_BATCHES_DIRNAME,
        root / "generic_index.sqlite",
    )
    if (
        normalized == "full"
        and not mode_dir.exists()
        and any(path.exists() for path in legacy_full_markers)
    ):
        return root
    return mode_dir


@dataclass(frozen=True)
class RecommendationArtifactProgress:
    """One progress update from the recommendation artifact workflow."""

    phase: str
    source_file_count: int
    shard_count: int
    row_count: int
    message: str


class RecommendationArtifactBuildCancelled(RuntimeError):
    """Raised after safely stopping an artifact build."""


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _artifact_entry(path: Path, root: Path) -> Dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "relative_path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
    }


def _selected_paths(dataset_path: ConversionDatasetInput) -> tuple[Path, ...]:
    if isinstance(dataset_path, (str, Path)):
        return (Path(dataset_path),)
    return tuple(Path(value) for value in dataset_path)


def _validate_paths(
    dataset_path: ConversionDatasetInput,
    output_dir: Path,
) -> None:
    destination = output_dir.resolve()
    for selected in _selected_paths(dataset_path):
        source = selected.resolve()
        if source.is_dir() and (
            destination == source or source in destination.parents
        ):
            raise ValueError(
                "Output folder must be outside the source folder(s) so "
                "generated CSV files are not discovered as source data."
            )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_batch_name(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "-", value.strip()).strip(".-")
    if not normalized:
        raise ValueError("Batch name must contain a letter or number")
    return normalized


def _automatic_batch_name(dataset_path: ConversionDatasetInput) -> str:
    selected = tuple(
        sorted(
            discover_conversion_datasets(dataset_path),
            key=lambda path: str(path.resolve()).casefold(),
        )
    )
    if not selected:
        selected = tuple(
            sorted(
                _selected_paths(dataset_path),
                key=lambda path: str(path.resolve()).casefold(),
            )
        )
    identity = json.dumps(
        [str(path.resolve()).casefold() for path in selected],
        ensure_ascii=False,
        separators=(",", ":"),
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:10]
    stem = _safe_batch_name(selected[0].stem)[:36] if selected else "inputs"
    return f"batch-{stem}-{digest}"


def _source_artifact_location(
    source: Path,
    library_dir: Path,
    *,
    source_sha256: str,
    shard_size: int,
    conversion_mode: str,
    use_rxnmapper: bool,
) -> tuple[str, Path]:
    """Return the immutable cache identity and folder for one source version."""
    resolved = source.resolve()
    path_digest = hashlib.sha256(
        str(resolved).casefold().encode("utf-8")
    ).hexdigest()[:10]
    identity = {
        "source_path": str(resolved).casefold(),
        "source_sha256": source_sha256,
        "shard_size": shard_size,
        "conversion_mode": conversion_mode.casefold(),
        "use_rxnmapper": use_rxnmapper,
        "sharded_conversion_definition_version": (
            SHARDED_CONVERSION_DEFINITION_VERSION
        ),
    }
    conversion_digest = hashlib.sha256(
        json.dumps(
            identity,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()[:16]
    try:
        stem = _safe_batch_name(source.stem)[:36]
    except ValueError:
        stem = "source"
    artifact_id = f"{stem}.{path_digest}/{conversion_digest}"
    return artifact_id, library_dir / CONVERTED_SOURCES_DIRNAME / artifact_id


def _relative_reference(path: Path, origin: Path) -> str:
    return Path(os.path.relpath(path.resolve(), origin.resolve())).as_posix()


def _resolve_manifest_reference(
    batch_manifest_path: Path,
    reference: Mapping[str, Any],
) -> Path:
    relative = str(reference.get("relative_path") or "").strip()
    if relative:
        return (batch_manifest_path.parent / relative).resolve()
    stored = str(reference.get("path") or "").strip()
    if not stored:
        raise ValueError(
            f"Saved batch has an empty source-manifest reference: "
            f"{batch_manifest_path}"
        )
    return Path(stored).resolve()


def _referenced_conversion_manifests(
    manifest_path: Path,
    manifest: Optional[Mapping[str, Any]] = None,
) -> tuple[Path, ...]:
    payload = (
        dict(manifest)
        if manifest is not None
        else json.loads(manifest_path.read_text(encoding="utf-8"))
    )
    artifact_type = payload.get("artifact_type")
    if artifact_type == "generic_sharded_conversion":
        return (manifest_path,)
    if artifact_type != SAVED_BATCH_MANIFEST_ARTIFACT_TYPE:
        raise ValueError(f"Unsupported saved batch manifest: {manifest_path}")
    if payload.get("schema_version") != SAVED_BATCH_WORKFLOW_SCHEMA_VERSION:
        raise ValueError(f"Unsupported saved batch manifest schema: {manifest_path}")
    references = tuple(payload.get("source_manifests") or ())
    if not references:
        raise ValueError(f"Saved batch has no converted-source references: {manifest_path}")
    return tuple(
        _resolve_manifest_reference(manifest_path, reference)
        for reference in references
    )


def _conversion_manifest_complete(manifest_path: Path) -> bool:
    if not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    source_files = tuple(manifest.get("source_files") or ())
    if (
        manifest.get("schema_version") != SHARD_MANIFEST_SCHEMA_VERSION
        or manifest.get("artifact_type") != "generic_sharded_conversion"
        or not source_files
        or any(not source.get("coverage_complete") for source in source_files)
        or any(
            entry.get("status") != "complete"
            for entry in manifest.get("shards") or ()
        )
    ):
        return False
    try:
        return bool(
            validate_sharded_conversion(manifest_path, verify_rows=False).get("valid")
        )
    except (OSError, ValueError, KeyError, TypeError):
        return False


def discover_saved_conversion_batches(
    library_dir: str | Path,
) -> tuple[Path, ...]:
    """Return root-level legacy and named saved-batch manifests."""
    root = Path(library_dir)
    manifests = []
    legacy = root / "shard_manifest.json"
    if legacy.is_file():
        manifests.append(legacy)
    batches_dir = root / SAVED_BATCHES_DIRNAME
    if batches_dir.is_dir():
        manifests.extend(
            sorted(
                batches_dir.glob("*/shard_manifest.json"),
                key=lambda path: path.parent.name.casefold(),
            )
        )
    return tuple(manifests)


def incomplete_saved_conversion_batches(
    library_dir: str | Path,
) -> tuple[Path, ...]:
    """Return saved manifests whose selected sources are not fully covered."""
    incomplete = []
    for manifest_path in discover_saved_conversion_batches(library_dir):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            incomplete.append(manifest_path)
            continue
        source_files = tuple(manifest.get("source_files") or ())
        if not source_files or any(
            not source.get("coverage_complete") for source in source_files
        ):
            incomplete.append(manifest_path)
            continue
        if manifest.get("artifact_type") == SAVED_BATCH_MANIFEST_ARTIFACT_TYPE:
            try:
                referenced = _referenced_conversion_manifests(
                    manifest_path, manifest
                )
            except ValueError:
                incomplete.append(manifest_path)
                continue
            if any(not _conversion_manifest_complete(path) for path in referenced):
                incomplete.append(manifest_path)
    return tuple(incomplete)


def resume_saved_conversion_batch(
    manifest_path: str | Path,
    *,
    workers: int = 1,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Resume one checkpointed batch using its persisted source selection."""
    source = Path(manifest_path)
    manifest = json.loads(source.read_text(encoding="utf-8"))
    artifact_type = manifest.get("artifact_type")
    if artifact_type == SAVED_BATCH_MANIFEST_ARTIFACT_TYPE:
        if manifest.get("schema_version") != SAVED_BATCH_WORKFLOW_SCHEMA_VERSION:
            raise ValueError(f"Unsupported saved batch manifest: {source}")
        dataset_paths = tuple(
            str(value) for value in manifest.get("dataset_paths") or ()
        )
        if not dataset_paths:
            raise ValueError(
                f"Saved batch has no source selection to resume: {source}"
            )
        missing = tuple(value for value in dataset_paths if not Path(value).exists())
        if missing:
            raise ValueError(
                f"Cannot resume {source}; {len(missing)} selected source path(s) "
                "no longer exist"
            )
        library_dir = Path(
            str(manifest.get("library_dir") or source.parent.parent)
        )
        return save_recommendation_batch(
            dataset_paths,
            library_dir,
            batch_name=str(manifest.get("batch_name") or source.parent.name),
            shard_size=int(manifest.get("shard_size") or 1_000),
            workers=workers,
            use_rxnmapper=bool(manifest.get("use_rxnmapper")),
            conversion_mode=str(manifest.get("mode") or "full"),
            progress_callback=progress_callback,
            cancel_check=cancel_check,
        )
    if (
        manifest.get("schema_version") != SHARD_MANIFEST_SCHEMA_VERSION
        or artifact_type != "generic_sharded_conversion"
    ):
        raise ValueError(f"Unsupported saved batch manifest: {source}")
    dataset_paths = tuple(str(value) for value in manifest.get("dataset_paths") or ())
    if not dataset_paths:
        raise ValueError(f"Saved batch has no source selection to resume: {source}")
    missing = tuple(value for value in dataset_paths if not Path(value).exists())
    if missing:
        raise ValueError(
            f"Cannot resume {source}; {len(missing)} selected source path(s) "
            "no longer exist"
        )
    mapping_contract = (
        (manifest.get("definition_contract") or {}).get("external_atom_mapping")
        or {}
    )
    return build_recommendation_artifacts(
        dataset_paths,
        source.parent,
        shard_size=int(manifest.get("shard_size") or 1_000),
        workers=workers,
        build_fast_index=False,
        use_rxnmapper=bool(mapping_contract.get("enabled")),
        conversion_mode=str(manifest.get("mode") or "full"),
        progress_callback=progress_callback,
        cancel_check=cancel_check,
    )


def _source_reference(
    *,
    source: Path,
    source_sha256: str,
    artifact_id: str,
    source_manifest: Path,
    batch_dir: Path,
) -> Dict[str, Any]:
    return {
        "source_path": str(source.resolve()),
        "source_sha256": source_sha256,
        "source_artifact_id": artifact_id,
        "path": str(source_manifest.resolve()),
        "relative_path": _relative_reference(source_manifest, batch_dir),
    }


def _saved_batch_manifest_payload(
    *,
    batch_name: str,
    batch_dir: Path,
    library_dir: Path,
    dataset_path: ConversionDatasetInput,
    source_references: Iterable[Mapping[str, Any]],
    shard_size: int,
    conversion_mode: str,
    use_rxnmapper: bool,
) -> Dict[str, Any]:
    references = [dict(reference) for reference in source_references]
    source_files = []
    for reference in references:
        manifest_path = _resolve_manifest_reference(
            batch_dir / "shard_manifest.json", reference
        )
        if not manifest_path.is_file():
            source_files.append(
                {
                    "path": reference["source_path"],
                    "sha256": reference["source_sha256"],
                    "covered_row_count": 0,
                    "coverage_complete": False,
                    "source_artifact_id": reference["source_artifact_id"],
                }
            )
            continue
        try:
            source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            source_files.append(
                {
                    "path": reference["source_path"],
                    "sha256": reference["source_sha256"],
                    "covered_row_count": 0,
                    "coverage_complete": False,
                    "source_artifact_id": reference["source_artifact_id"],
                }
            )
            continue
        for source_file in source_manifest.get("source_files") or ():
            source_files.append(
                {
                    **source_file,
                    "source_artifact_id": reference["source_artifact_id"],
                }
            )
    selected = _selected_paths(dataset_path)
    return {
        "schema_version": SAVED_BATCH_WORKFLOW_SCHEMA_VERSION,
        "artifact_type": SAVED_BATCH_MANIFEST_ARTIFACT_TYPE,
        "batch_name": batch_name,
        "batch_dir": str(batch_dir.resolve()),
        "library_dir": str(library_dir.resolve()),
        "dataset_path": (
            str(selected[0].resolve()) if len(selected) == 1 else None
        ),
        "dataset_paths": [str(path.resolve()) for path in selected],
        "mode": conversion_mode,
        "shard_size": shard_size,
        "use_rxnmapper": use_rxnmapper,
        "source_files": source_files,
        "source_manifests": references,
        "coverage_complete": bool(source_files)
        and all(source.get("coverage_complete") for source in source_files),
    }


def save_recommendation_batch(
    dataset_path: ConversionDatasetInput,
    library_dir: str | Path,
    *,
    batch_name: str = "",
    shard_size: int = 1_000,
    workers: int = 1,
    use_rxnmapper: bool = False,
    conversion_mode: str = "full",
    checkpoint_interval: int = 1,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Convert sources once and save a lightweight batch of artifact references."""
    root = Path(library_dir)
    batches_dir = root / SAVED_BATCHES_DIRNAME
    resolved_name = (
        _safe_batch_name(batch_name)
        if batch_name.strip()
        else _automatic_batch_name(dataset_path)
    )
    batch_dir = batches_dir / resolved_name
    _validate_paths(dataset_path, batch_dir)
    batches_dir.mkdir(parents=True, exist_ok=True)
    batch_dir.mkdir(parents=True, exist_ok=True)
    paths = discover_conversion_datasets(dataset_path)
    if not paths:
        raise ValueError("No supported conversion datasets were found")

    source_references = []
    source_work = []
    for path in paths:
        source_sha256 = _sha256(path)
        artifact_id, artifact_dir = _source_artifact_location(
            path,
            root,
            source_sha256=source_sha256,
            shard_size=shard_size,
            conversion_mode=conversion_mode,
            use_rxnmapper=use_rxnmapper,
        )
        source_manifest = artifact_dir / "shard_manifest.json"
        source_references.append(
            _source_reference(
                source=path,
                source_sha256=source_sha256,
                artifact_id=artifact_id,
                source_manifest=source_manifest,
                batch_dir=batch_dir,
            )
        )
        source_work.append((path, artifact_dir, source_manifest))

    batch_manifest_path = batch_dir / "shard_manifest.json"

    def checkpoint_batch() -> Dict[str, Any]:
        payload = _saved_batch_manifest_payload(
            batch_name=resolved_name,
            batch_dir=batch_dir,
            library_dir=root,
            dataset_path=dataset_path,
            source_references=source_references,
            shard_size=shard_size,
            conversion_mode=conversion_mode,
            use_rxnmapper=use_rxnmapper,
        )
        _atomic_json(batch_manifest_path, payload)
        return payload

    checkpoint_batch()
    source_reports = []
    source_reused = []
    reused_source_file_count = 0
    try:
        for path, artifact_dir, source_manifest in source_work:
            if _conversion_manifest_complete(source_manifest):
                stored_report_path = artifact_dir / "recommendation_artifacts_report.json"
                if stored_report_path.is_file():
                    source_report = json.loads(
                        stored_report_path.read_text(encoding="utf-8")
                    )
                    source_report = {
                        **source_report,
                        "reused_shard_count": int(source_report.get("shard_count") or 0),
                    }
                    reused_source_file_count += 1
                    source_reports.append(source_report)
                    source_reused.append(True)
                    if progress_callback is not None:
                        progress_callback(
                            RecommendationArtifactProgress(
                                phase="source_reused",
                                source_file_count=len(paths),
                                shard_count=int(
                                    source_report.get("shard_count") or 0
                                ),
                                row_count=int(
                                    source_report.get("record_count") or 0
                                ),
                                message=(
                                    f"Reused converted source: {path.name}"
                                ),
                            )
                        )
                    checkpoint_batch()
                    continue
            source_report = build_recommendation_artifacts(
                path,
                artifact_dir,
                shard_size=shard_size,
                workers=workers,
                build_fast_index=False,
                use_rxnmapper=use_rxnmapper,
                conversion_mode=conversion_mode,
                checkpoint_interval=checkpoint_interval,
                progress_callback=progress_callback,
                cancel_check=cancel_check,
            )
            source_reports.append(source_report)
            source_reused.append(False)
            checkpoint_batch()
    except RecommendationArtifactBuildCancelled:
        checkpoint_batch()
        raise

    batch_manifest = checkpoint_batch()
    record_count = sum(int(item.get("record_count") or 0) for item in source_reports)
    shard_count = sum(int(item.get("shard_count") or 0) for item in source_reports)
    reused_shard_count = sum(
        int(item.get("reused_shard_count") or 0) for item in source_reports
    )
    source_artifacts = [
        {
            **reference,
            "record_count": int(report.get("record_count") or 0),
            "shard_count": int(report.get("shard_count") or 0),
            "reused": reused,
        }
        for reference, report, reused in zip(
            source_references, source_reports, source_reused
        )
    ]
    report = {
        "schema_version": RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION,
        "artifact_type": "recommendation_artifact_build",
        "conversion_mode": conversion_mode,
        "source_path": (
            str(paths[0].resolve()) if len(paths) == 1 else None
        ),
        "source_paths": [str(path.resolve()) for path in paths],
        "output_dir": str(batch_dir.resolve()),
        "settings": {
            "shard_size": shard_size,
            "workers": 1 if use_rxnmapper else workers,
            "requested_workers": workers,
            "build_fast_index": False,
            "use_rxnmapper": use_rxnmapper,
            "compression": "shared gzip canonical shards",
        },
        "source_file_count": len(paths),
        "record_count": record_count,
        "shard_count": shard_count,
        "reused_shard_count": reused_shard_count,
        "reused_source_file_count": reused_source_file_count,
        "new_source_file_count": len(paths) - reused_source_file_count,
        "artifacts": {
            "canonical_manifest": _artifact_entry(batch_manifest_path, batch_dir),
            "converted_sources": source_artifacts,
        },
        "warnings": [],
        "batch_manifest": batch_manifest,
    }
    batch_report = {
        **report,
        "artifact_type": "saved_recommendation_batch",
        "batch_name": resolved_name,
        "batch_dir": str(batch_dir.resolve()),
        "library_dir": str(root.resolve()),
        "conversion_mode": conversion_mode,
    }
    _atomic_json(batch_dir / "saved_batch_report.json", batch_report)
    return batch_report


def _batch_records(
    manifest_paths: Iterable[Path],
) -> Iterator[Dict[str, Any]]:
    expected_contract: Mapping[str, Any] | None = None
    for batch_manifest_path in manifest_paths:
        batch_manifest = json.loads(
            batch_manifest_path.read_text(encoding="utf-8")
        )
        source_files = tuple(batch_manifest.get("source_files") or ())
        incomplete_sources = tuple(
            str(source.get("path") or "")
            for source in source_files
            if not source.get("coverage_complete")
        )
        if not source_files or incomplete_sources:
            raise ValueError(
                "Saved batch conversion is incomplete and cannot be combined: "
                f"{batch_manifest_path} ({len(incomplete_sources)} incomplete "
                "source file(s)). Resume and finish this batch first."
            )
        conversion_manifests = _referenced_conversion_manifests(
            batch_manifest_path, batch_manifest
        )
        for manifest_path in conversion_manifests:
            if not manifest_path.is_file():
                raise ValueError(
                    "Saved batch references a missing converted source: "
                    f"{manifest_path}"
                )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if (
                manifest.get("schema_version") != SHARD_MANIFEST_SCHEMA_VERSION
                or manifest.get("artifact_type") != "generic_sharded_conversion"
            ):
                raise ValueError(
                    f"Unsupported converted-source manifest: {manifest_path}"
                )
            contract = {
                key: value
                for key, value in (manifest.get("definition_contract") or {}).items()
                if key != "external_atom_mapping"
            }
            if expected_contract is None:
                expected_contract = contract
            elif contract != expected_contract:
                raise ValueError(
                    "Saved batches use different chemistry definition contracts; "
                    "reconvert the older batch before combining"
                )
            for entry in manifest.get("shards") or ():
                if entry.get("status") != "complete":
                    raise ValueError(
                        "Saved batch contains an incomplete shard: "
                        f"{manifest_path}"
                    )
                shard_path = manifest_path.parent / str(
                    entry.get("output_path") or ""
                )
                if not shard_path.is_file() or _sha256(shard_path) != entry.get(
                    "output_sha256"
                ):
                    raise ValueError(
                        f"Saved batch shard checksum failed: {shard_path}"
                    )
                yield from iter_gzip_jsonl(shard_path)


def _write_combined_records(
    manifest_paths: tuple[Path, ...],
    output_path: Path,
    *,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    dedupe_path = output_path.with_suffix(output_path.suffix + ".dedupe.tmp.sqlite")
    input_count = 0
    output_count = 0
    duplicate_count = 0
    precedent_tier_counts: Dict[str, int] = {}
    core_eligibility_counts: Dict[str, int] = {}
    if dedupe_path.is_file():
        dedupe_path.unlink()
    dedupe = sqlite3.connect(dedupe_path)
    try:
        dedupe.execute(
            "CREATE TABLE seen_records ("
            "identity TEXT PRIMARY KEY, digest TEXT NOT NULL) WITHOUT ROWID"
        )
        with temporary.open("wb") as raw:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw,
                compresslevel=6,
                mtime=0,
            ) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8") as text:
                    for record in _batch_records(manifest_paths):
                        if cancel_check is not None and cancel_check():
                            raise RecommendationArtifactBuildCancelled(
                                "Combine cancelled while writing canonical records"
                            )
                        input_count += 1
                        canonical = json.dumps(
                            record,
                            ensure_ascii=False,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
                        observation_id = str(record.get("observation_id") or "")
                        identity = f"observation:{observation_id}" if observation_id else (
                            f"record:{digest}"
                        )
                        previous_row = dedupe.execute(
                            "SELECT digest FROM seen_records WHERE identity = ?",
                            (identity,),
                        ).fetchone()
                        if previous_row is not None:
                            previous = str(previous_row[0])
                            if previous != digest:
                                raise ValueError(
                                    "Conflicting saved records share observation ID "
                                    f"{observation_id!r}"
                                )
                            duplicate_count += 1
                            continue
                        dedupe.execute(
                            "INSERT INTO seen_records VALUES (?, ?)",
                            (identity, digest),
                        )
                        text.write(canonical)
                        text.write("\n")
                        output_count += 1
                        precedent = str(record.get("precedent_tier") or "unavailable")
                        precedent_tier_counts[precedent] = (
                            precedent_tier_counts.get(precedent, 0) + 1
                        )
                        core = str(record.get("core_eligibility") or "")
                        core_eligibility_counts[core] = (
                            core_eligibility_counts.get(core, 0) + 1
                        )
                        if output_count % 1_000 == 0:
                            dedupe.commit()
                        if output_count % 1_000 == 0 and progress_callback is not None:
                            progress_callback(
                                RecommendationArtifactProgress(
                                    phase="combine_records",
                                    source_file_count=len(manifest_paths),
                                    shard_count=0,
                                    row_count=output_count,
                                    message=(
                                        f"Combined {output_count} unique reaction(s)."
                                    ),
                                )
                            )
        temporary.replace(output_path)
    finally:
        dedupe.close()
        if temporary.is_file():
            temporary.unlink()
        if dedupe_path.is_file():
            dedupe_path.unlink()
    return {
        "input_record_count": input_count,
        "record_count": output_count,
        "duplicate_record_count": duplicate_count,
        "precedent_tier_counts": dict(sorted(precedent_tier_counts.items())),
        "core_eligibility_counts": dict(sorted(core_eligibility_counts.items())),
    }


def combine_saved_recommendation_batches(
    library_dir: str | Path,
    *,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Combine every saved batch and build the active recommender artifacts."""
    destination = Path(library_dir)
    destination.mkdir(parents=True, exist_ok=True)
    manifests = discover_saved_conversion_batches(destination)
    if not manifests:
        raise ValueError("No saved conversion batches were found")

    def notify(phase: str, message: str, row_count: int = 0) -> None:
        if progress_callback is not None:
            progress_callback(
                RecommendationArtifactProgress(
                    phase=phase,
                    source_file_count=len(manifests),
                    shard_count=0,
                    row_count=row_count,
                    message=message,
                )
            )

    notify("combine_started", f"Combining {len(manifests)} saved batch(es)…")
    records_path = destination / COMBINED_RECORDS_FILENAME
    combined = _write_combined_records(
        manifests,
        records_path,
        progress_callback=progress_callback,
        cancel_check=cancel_check,
    )
    total_rows = int(combined["record_count"])
    if cancel_check is not None and cancel_check():
        raise RecommendationArtifactBuildCancelled(
            "Combine cancelled before catalog creation"
        )
    notify("catalogs_started", "Building combined condition and reference catalogs…")
    catalogs = write_conversion_catalogs(
        iter_canonical_records(records_path), destination
    )

    index_path = destination / "generic_index.sqlite"
    review_index_path = destination / "generic_review_index.sqlite"

    def on_index_progress(phase: str, count: int) -> None:
        message = {
            "staged": f"Staged {count} eligible precedent(s) on disk.",
            "rows": f"Wrote {count} sorted precedent row(s) to the index.",
        }.get(phase, f"Aggregated {count} index lookup key(s).")
        notify("index_rows", message, count)

    notify("index_started", "Building the combined fast recommendation index…")
    try:
        index_report = build_sqlite_generic_index(
            iter_canonical_records(records_path),
            index_path,
            progress_callback=on_index_progress,
            cancel_check=cancel_check,
        )
        if review_index_path.is_file():
            review_index_path.unlink()
    except SQLiteIndexBuildCancelled as exc:
        raise RecommendationArtifactBuildCancelled(str(exc)) from exc

    combined_manifest_path = destination / COMBINED_BATCH_MANIFEST_FILENAME
    combined_manifest = {
        "schema_version": SAVED_BATCH_WORKFLOW_SCHEMA_VERSION,
        "artifact_type": "combined_saved_recommendation_batches",
        "batch_manifests": [
            {
                "path": str(path.resolve()),
                "sha256": _sha256(path),
            }
            for path in manifests
        ],
        "combined_records": {
            "path": records_path.name,
            "sha256": _sha256(records_path),
            "row_count": total_rows,
        },
        **combined,
    }
    _atomic_json(combined_manifest_path, combined_manifest)
    trusted_count = int(index_report["row_count"])
    review_core_count = int(
        combined["precedent_tier_counts"].get("review_core") or 0
    )
    core_counts = combined["core_eligibility_counts"]
    query_core_count = sum(
        int(core_counts.get(tier) or 0)
        for tier in ("trusted_core", "review_core", "query_only")
    )
    report: Dict[str, Any] = {
        "schema_version": SAVED_BATCH_WORKFLOW_SCHEMA_VERSION,
        "artifact_type": "combined_recommendation_artifact_build",
        "output_dir": str(destination.resolve()),
        "batch_count": len(manifests),
        "batch_manifest_paths": [str(path.resolve()) for path in manifests],
        **combined,
        "eligible_index_record_count": trusted_count,
        "trusted_precedent_count": trusted_count,
        "review_core_precedent_count": review_core_count,
        "unrestricted_precedent_count": None,
        "review_index_generated": False,
        "review_index_reuses_trusted": False,
        "query_core_eligible_count": query_core_count,
        "artifacts": {
            "combined_manifest": _artifact_entry(
                combined_manifest_path, destination
            ),
            "canonical_records": _artifact_entry(records_path, destination),
            "fast_index": _artifact_entry(index_path, destination),
        },
        "catalogs": catalogs,
        "warnings": [],
    }
    report_path = destination / "combined_recommendation_report.json"
    report["report_path"] = str(report_path.resolve())
    _atomic_json(report_path, report)
    notify(
        "completed",
        f"Combined index ready for {total_rows} unique reaction(s).",
        total_rows,
    )
    return report


def build_recommendation_artifacts(
    dataset_path: ConversionDatasetInput,
    output_dir: str | Path,
    *,
    shard_size: int = 1_000,
    workers: int = 1,
    build_fast_index: bool = True,
    use_rxnmapper: bool = False,
    conversion_mode: str = "full",
    checkpoint_interval: int = 1,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Build restartable canonical data and an optional fast index.

    Canonical records remain the source of truth. The persisted generic index
    contains only retrieval-eligible fields and lookup maps. Human review CSVs
    remain available through the standalone review-export API and CLI, but are
    not generated by the converter workflow.
    """
    destination = Path(output_dir)
    selected_paths = _selected_paths(dataset_path)
    _validate_paths(selected_paths, destination)
    destination.mkdir(parents=True, exist_ok=True)

    latest_source_file_count = 0
    latest_shard_count = 0

    def notify(
        phase: str,
        message: str,
        *,
        row_count: int = 0,
        source_file_count: Optional[int] = None,
        shard_count: Optional[int] = None,
    ) -> None:
        nonlocal latest_source_file_count, latest_shard_count
        if source_file_count is not None:
            latest_source_file_count = source_file_count
        if shard_count is not None:
            latest_shard_count = shard_count
        if progress_callback is not None:
            progress_callback(
                RecommendationArtifactProgress(
                    phase=phase,
                    source_file_count=latest_source_file_count,
                    shard_count=latest_shard_count,
                    row_count=row_count,
                    message=message,
                )
            )

    def on_conversion_progress(progress: ShardedConversionProgress) -> None:
        notify(
            f"canonical_{progress.phase}",
            progress.message,
            row_count=progress.row_count,
            source_file_count=progress.source_file_count,
            shard_count=progress.shard_count,
        )

    effective_workers = 1 if use_rxnmapper else workers
    try:
        conversion_report = convert_datasets_sharded(
            selected_paths,
            destination,
            shard_size=shard_size,
            mode=conversion_mode,
            workers=effective_workers,
            checkpoint_interval=checkpoint_interval,
            merge_records=False,
            use_rxnmapper=use_rxnmapper,
            progress_callback=on_conversion_progress,
            cancel_check=cancel_check,
        )
    except ShardedConversionCancelled as exc:
        raise RecommendationArtifactBuildCancelled(str(exc)) from exc

    failed_shards = int(conversion_report.get("failed_shard_count") or 0)
    if failed_shards:
        failure_entries = conversion_report.get("failed_shards") or ()
        first_failure = (
            (failure_entries[0].get("failures") or [{}])[0]
            if failure_entries
            else {}
        )
        first_detail = str(first_failure.get("message") or "").strip()
        detail_suffix = (
            f" First failure: {first_detail}" if first_detail else ""
        )
        raise RuntimeError(
            f"Canonical conversion has {failed_shards} failed shard(s); "
            "index artifacts were not built."
            f"{detail_suffix} See {destination / 'shard_manifest.json'}."
        )
    records_path = destination / "shard_manifest.json"
    total_rows = int(conversion_report.get("output_row_count") or 0)
    if cancel_check is not None and cancel_check():
        raise RecommendationArtifactBuildCancelled(
            "Build cancelled after canonical data completed; run again with "
            "the same output folder to reuse it."
        )

    index_report: Optional[Dict[str, Any]] = None
    index_path = destination / "generic_index.sqlite"
    review_index_path = destination / "generic_review_index.sqlite"
    retired_json_index_paths = (
        destination / "generic_index.json.gz",
        destination / "generic_review_index.json.gz",
    )
    for retired_path in retired_json_index_paths:
        if retired_path.is_file():
            retired_path.unlink()
    if build_fast_index:
        if cancel_check is not None and cancel_check():
            raise RecommendationArtifactBuildCancelled(
                "Build cancelled before index creation; canonical data are "
                "complete."
            )
        notify(
            "index_started",
            "Building the compressed fast-load recommendation index…",
            row_count=0,
        )
        scanned_rows = 0

        def records() -> Iterator[Dict[str, Any]]:
            nonlocal scanned_rows
            for record in iter_canonical_records(records_path):
                if cancel_check is not None and cancel_check():
                    raise RecommendationArtifactBuildCancelled(
                        "Build cancelled during index creation; canonical data "
                        "are complete."
                    )
                scanned_rows += 1
                if scanned_rows % 1_000 == 0:
                    notify(
                        "index_rows",
                        f"Indexed input scan: {scanned_rows} record(s).",
                        row_count=scanned_rows,
                    )
                yield record

        def on_index_progress(phase: str, count: int) -> None:
            if phase == "staged":
                message = f"Staged {count} eligible precedent(s) on disk."
            elif phase == "rows":
                message = f"Wrote {count} sorted precedent row(s) to the index."
            else:
                message = f"Aggregated {count} index lookup key(s)."
            notify("index_rows", message, row_count=count)

        try:
            index_report = build_sqlite_generic_index(
                records(),
                index_path,
                progress_callback=on_index_progress,
                cancel_check=cancel_check,
            )
            if review_index_path.is_file():
                review_index_path.unlink()
        except SQLiteIndexBuildCancelled as exc:
            raise RecommendationArtifactBuildCancelled(str(exc)) from exc
        notify(
            "index_completed",
            (
                "Fast-load index complete: "
                f"{index_report['row_count']} trusted precedent(s)."
            ),
            row_count=int(index_report["row_count"]),
        )

    artifacts: Dict[str, Any] = {
        "canonical_manifest": _artifact_entry(records_path, destination),
    }
    if index_report is not None:
        artifacts["fast_index"] = _artifact_entry(index_path, destination)
    shard_paths = tuple((destination / "shards").glob("*.jsonl.gz"))
    shard_size_bytes = sum(path.stat().st_size for path in shard_paths)
    all_files = tuple(path for path in destination.rglob("*") if path.is_file())
    warnings = []
    if use_rxnmapper and workers != effective_workers:
        warnings.append(
            "RXNMapper uses one conversion worker to avoid loading a separate "
            "model in each process."
        )
    if not build_fast_index and any(
        path.is_file()
        for path in (index_path, review_index_path)
    ):
        warnings.append(
            "Older trusted or review-core indexes may exist but were not rebuilt; "
            "do not use them unless they match the current canonical records."
        )
    core_eligibility_counts = dict(
        conversion_report.get("core_eligibility_counts") or {}
    )
    trusted_precedent_count = (
        int(index_report["row_count"]) if index_report is not None else None
    )
    review_core_count = int(
        (conversion_report.get("precedent_tier_counts") or {}).get(
            "review_core", 0
        )
    )
    query_core_count = sum(
        int(core_eligibility_counts.get(tier) or 0)
        for tier in ("trusted_core", "review_core", "query_only")
    )
    report: Dict[str, Any] = {
        "schema_version": RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION,
        "artifact_type": "recommendation_artifact_build",
        "conversion_mode": conversion_mode,
        "source_path": (
            str(selected_paths[0].resolve()) if len(selected_paths) == 1 else None
        ),
        "source_paths": [str(path.resolve()) for path in selected_paths],
        "output_dir": str(destination.resolve()),
        "settings": {
            "shard_size": shard_size,
            "workers": effective_workers,
            "requested_workers": workers,
            "build_fast_index": build_fast_index,
            "use_rxnmapper": use_rxnmapper,
            "compression": "gzip canonical shards; row-compressed SQLite index",
        },
        "source_file_count": int(conversion_report["source_file_count"]),
        "record_count": total_rows,
        "eligible_index_record_count": (
            trusted_precedent_count
        ),
        "trusted_precedent_count": trusted_precedent_count,
        "review_core_precedent_count": review_core_count,
        "review_index_generated": False,
        "review_index_reuses_trusted": False,
        "unrestricted_precedent_count": None,
        "query_core_eligible_count": query_core_count,
        "core_eligibility_counts": core_eligibility_counts,
        "shard_count": int(conversion_report["shard_count"]),
        "reused_shard_count": int(conversion_report["reused_shard_count"]),
        "artifacts": artifacts,
        "storage": {
            "shard_file_count": len(shard_paths),
            "shard_size_bytes": shard_size_bytes,
            "total_output_size_bytes_before_report": sum(
                path.stat().st_size for path in all_files
            ),
        },
        "warnings": warnings,
        "conversion_report_path": str(
            (destination / "conversion_report.json").resolve()
        ),
    }
    report_path = destination / "recommendation_artifacts_report.json"
    report["report_path"] = str(report_path.resolve())
    _atomic_json(report_path, report)
    notify(
        "completed",
        (
            f"All requested artifacts are ready for {total_rows} reaction(s)."
        ),
        row_count=total_rows,
    )
    return report


__all__ = [
    "COMBINED_BATCH_MANIFEST_FILENAME",
    "COMBINED_RECORDS_FILENAME",
    "CONVERTED_SOURCES_DIRNAME",
    "RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION",
    "SAVED_BATCHES_DIRNAME",
    "SAVED_BATCH_MANIFEST_ARTIFACT_TYPE",
    "SAVED_BATCH_WORKFLOW_SCHEMA_VERSION",
    "RecommendationArtifactBuildCancelled",
    "RecommendationArtifactProgress",
    "build_recommendation_artifacts",
    "combine_saved_recommendation_batches",
    "discover_saved_conversion_batches",
    "incomplete_saved_conversion_batches",
    "resume_saved_conversion_batch",
    "save_recommendation_batch",
]
