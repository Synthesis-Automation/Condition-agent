"""Restartable, content-addressed conversion for the full reaction corpus."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import shutil
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from enum import Enum
from itertools import islice
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Iterator,
    Mapping,
    Optional,
    Sequence,
)

from condition_registry import condition_registry_definition_versions
from reactive_taxonomy import (
    REACTION_SIGNATURE_SCHEMA_VERSION,
    AtomMappingProvider,
    RxnMapperProvider,
    reaction_signature_definition_versions,
)
from rdkit import RDLogger

from ..models import (
    GENERIC_CONVERTER_DEFINITION_VERSION,
    RECOMMENDATION_RECORD_SCHEMA_VERSION,
)
from .generic import GenericConversionCache, convert_record
from .input_schema import RawReactionRecord, discover_csv_datasets, iter_csv_records

SHARD_MANIFEST_SCHEMA_VERSION = "1.0"
SHARDED_CONVERSION_DEFINITION_VERSION = "generic_sharded_conversion.v4.0"


@dataclass(frozen=True)
class ShardedConversionProgress:
    """One progress update from restartable sharded conversion."""

    phase: str
    source_file_count: int
    shard_count: int
    row_count: int
    message: str


class ShardedConversionCancelled(RuntimeError):
    """Raised after completed shards are checkpointed for a cancelled run."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _definition_contract(
    mapping_provider: AtomMappingProvider | None = None,
) -> Dict[str, Any]:
    mapping_metadata = (
        mapping_provider.metadata if mapping_provider is not None else None
    )
    return {
        "record_schema_version": RECOMMENDATION_RECORD_SCHEMA_VERSION,
        "converter_definition_version": GENERIC_CONVERTER_DEFINITION_VERSION,
        "reaction_signature_schema_version": REACTION_SIGNATURE_SCHEMA_VERSION,
        "taxonomy_definition_versions": reaction_signature_definition_versions(),
        "condition_registry_definition_versions": (
            condition_registry_definition_versions()
        ),
        "sharded_conversion_definition_version": (
            SHARDED_CONVERSION_DEFINITION_VERSION
        ),
        "external_atom_mapping": {
            "enabled": mapping_provider is not None,
            "provider_id": (
                mapping_metadata.provider_id if mapping_metadata is not None else None
            ),
            "provider_version": (
                mapping_metadata.provider_version
                if mapping_metadata is not None
                else None
            ),
            "model_id": (
                mapping_metadata.model_id if mapping_metadata is not None else None
            ),
            "model_sha256": (
                mapping_metadata.model_sha256
                if mapping_metadata is not None
                else None
            ),
        },
    }


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _write_gzip_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> int:
    temporary = path.with_suffix(path.suffix + ".tmp")
    count = 0
    with temporary.open("wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw,
            compresslevel=6,
            mtime=0,
        ) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8") as text:
                for row in rows:
                    text.write(
                        json.dumps(
                            row,
                            ensure_ascii=False,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                    )
                    text.write("\n")
                    count += 1
    temporary.replace(path)
    return count


def iter_gzip_jsonl(path: str | Path) -> Iterator[Dict[str, Any]]:
    """Stream canonical objects from one deterministic compressed shard."""
    with gzip.open(Path(path), mode="rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid shard JSONL at {path}:{line_number}: {exc.msg}"
                ) from exc
            if not isinstance(value, dict):
                raise ValueError(f"Shard row is not an object at {path}:{line_number}")
            yield value


def _chunks(
    records: Iterable[RawReactionRecord],
    shard_size: int,
) -> Iterator[tuple[RawReactionRecord, ...]]:
    iterator = iter(records)
    while True:
        values = tuple(islice(iterator, shard_size))
        if not values:
            return
        yield values


def _source_key(path: Path) -> str:
    digest = hashlib.sha256(str(path.resolve()).encode("utf-8")).hexdigest()[:10]
    safe_stem = "".join(
        character if character.isalnum() or character in {"-", "_"} else "_"
        for character in path.stem
    )
    return f"{safe_stem}.{digest}"


def _shard_id(path: Path, part_number: int) -> str:
    return f"{_source_key(path)}.part-{part_number:05d}"


def _reusable(
    previous: Mapping[str, Any],
    expected: Mapping[str, Any],
    destination: Path,
) -> bool:
    if previous.get("status") != "complete":
        return False
    keys = (
        "source_path",
        "source_sha256",
        "part_number",
        "input_position_start",
        "input_position_end",
        "input_row_count",
        "definition_contract",
        "output_path",
    )
    if any(previous.get(key) != expected.get(key) for key in keys):
        return False
    output = destination / str(previous["output_path"])
    return output.is_file() and _sha256(output) == previous.get("output_sha256")


def _converted_counts(payloads: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    fields = (
        "admission_tier",
        "chemistry_status",
        "condition_status",
        "condition_stage_status",
        "outcome_status",
        "index_eligibility",
        "evidence_quality",
    )

    def text(value: Any) -> str:
        return str(value.value if isinstance(value, Enum) else value or "")

    values = {
        field: Counter(text(payload.get(field)) for payload in payloads)
        for field in fields
    }
    return {
        f"{field}_counts": dict(sorted(counts.items()))
        for field, counts in values.items()
    } | {
        "signature_count": sum(
            int(payload.get("reaction_signature") is not None)
            for payload in payloads
        ),
        "admission_reason_counts": dict(
            sorted(
                Counter(
                    str(reason)
                    for payload in payloads
                    for reason in payload.get("admission_reasons") or ()
                ).items()
            )
        ),
        "transformation_class_counts": dict(
            sorted(
                Counter(
                    str(payload.get("transformation_class") or "unknown")
                    for payload in payloads
                ).items()
            )
        ),
        "named_family_counts": dict(
            sorted(
                Counter(
                    str(payload.get("named_family") or "unnamed")
                    for payload in payloads
                ).items()
            )
        ),
        "reaction_scope_counts": dict(
            sorted(
                Counter(
                    str(
                        (
                            (payload.get("reaction_signature") or {}).get(
                                "topology"
                            )
                            or {}
                        ).get("reaction_scope")
                        or "unknown"
                    )
                    for payload in payloads
                ).items()
            )
        ),
        "reaction_completeness_status_counts": dict(
            sorted(
                Counter(
                    str(
                        (
                            payload.get("reaction_completeness") or {}
                        ).get("status")
                        or "missing"
                    )
                    for payload in payloads
                ).items()
            )
        ),
        "external_mapping_status_counts": dict(
            sorted(
                Counter(
                    str(mapping.get("status") or "missing")
                    for payload in payloads
                    for mapping in (payload.get("external_atom_mapping"),)
                    if isinstance(mapping, Mapping)
                ).items()
            )
        ),
        "fragment_source_support_status_counts": dict(
            sorted(
                Counter(
                    str(support.get("status") or "missing")
                    for payload in payloads
                    for support in payload.get("fragment_source_support") or ()
                    if isinstance(support, Mapping)
                ).items()
            )
        ),
    }


def _convert_shard_task(
    task: tuple[Mapping[str, Any], tuple[RawReactionRecord, ...], str],
    *,
    cache: GenericConversionCache | None = None,
    mapping_provider: AtomMappingProvider | None = None,
) -> Dict[str, Any]:
    expected, raw_records, output_value = task
    output = Path(output_value)
    RDLogger.DisableLog("rdApp.warning")
    RDLogger.DisableLog("rdApp.error")
    conversion_cache = cache or GenericConversionCache(
        max_entries=max(1_000, len(raw_records) * 2)
    )
    try:
        payloads = [
            convert_record(
                raw_record,
                cache=conversion_cache,
                mapping_provider=mapping_provider,
            ).to_dict()
            for raw_record in raw_records
        ]
        output_count = _write_gzip_jsonl(output, payloads)
        return {
            **expected,
            "status": "complete",
            "output_row_count": output_count,
            "output_sha256": _sha256(output),
            "failure_count": 0,
            "failures": [],
            "reused": False,
            **_converted_counts(payloads),
        }
    except Exception as exc:
        return {
            **expected,
            "status": "failed",
            "output_row_count": 0,
            "output_sha256": "",
            "failure_count": 1,
            "failures": [
                {
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                }
            ],
            "reused": False,
        }


def _manifest_payload(
    *,
    dataset_path: str | Path,
    mode: str,
    shard_size: int,
    max_shards: Optional[int],
    contract: Mapping[str, Any],
    paths: Sequence[Path],
    source_checksums: Mapping[str, str],
    source_row_counts: Mapping[str, int],
    entries: Sequence[Mapping[str, Any]],
    coverage_complete: bool,
) -> Dict[str, Any]:
    return {
        "schema_version": SHARD_MANIFEST_SCHEMA_VERSION,
        "artifact_type": "generic_sharded_conversion",
        "dataset_path": str(Path(dataset_path).resolve()),
        "mode": mode,
        "shard_size": shard_size,
        "max_shards": max_shards,
        "definition_contract": contract,
        "source_files": [
            {
                "path": str(path.resolve()),
                "sha256": source_checksums[str(path.resolve())],
                "covered_row_count": source_row_counts.get(
                    str(path.resolve()), 0
                ),
                "coverage_complete": coverage_complete,
            }
            for path in paths
        ],
        "shards": sorted(
            entries,
            key=lambda entry: (
                str(entry["source_path"]).casefold(),
                int(entry["part_number"]),
            ),
        ),
    }
def _merge_counts(
    entries: Iterable[Mapping[str, Any]],
    field: str,
) -> Dict[str, int]:
    counts = Counter()
    for entry in entries:
        counts.update(entry.get(field) or {})
    return dict(sorted(counts.items()))


def _write_catalogs(
    manifest: Mapping[str, Any],
    destination: Path,
) -> Dict[str, Any]:
    recipes: Dict[str, Mapping[str, Any]] = {}
    references: Dict[str, Mapping[str, Any]] = {}
    series: Dict[str, Dict[str, Any]] = {}
    for entry in manifest["shards"]:
        for record in iter_gzip_jsonl(destination / entry["output_path"]):
            recipe_id = str(record.get("resolved_recipe_id") or "")
            if recipe_id and isinstance(record.get("resolved_recipe"), Mapping):
                recipes.setdefault(recipe_id, record["resolved_recipe"])
            reference_id = str(record.get("reference_id") or "")
            if reference_id and isinstance(record.get("reference_identity"), Mapping):
                references.setdefault(reference_id, record["reference_identity"])
            series_id = str(record.get("reference_condition_series_id") or "")
            if series_id:
                item = series.setdefault(
                    series_id,
                    {
                        "reference_condition_series_id": series_id,
                        "reference_id": reference_id,
                        "recipe_core_ids": set(),
                        "canonical_reaction_ids": set(),
                        "observation_count": 0,
                    },
                )
                recipe_core_id = str(record.get("resolved_recipe_core_id") or "")
                reaction_id = str(record.get("canonical_reaction_id") or "")
                if recipe_core_id:
                    item["recipe_core_ids"].add(recipe_core_id)
                if reaction_id:
                    item["canonical_reaction_ids"].add(reaction_id)
                item["observation_count"] += 1
    paths = {
        "recipe_catalog": Path("recipe_catalog.jsonl.gz"),
        "reference_catalog": Path("reference_catalog.jsonl.gz"),
        "reference_condition_series_catalog": Path(
            "reference_condition_series.jsonl.gz"
        ),
    }
    _write_gzip_jsonl(
        destination / paths["recipe_catalog"],
        (recipes[key] for key in sorted(recipes)),
    )
    _write_gzip_jsonl(
        destination / paths["reference_catalog"],
        (references[key] for key in sorted(references)),
    )
    _write_gzip_jsonl(
        destination / paths["reference_condition_series_catalog"],
        (
            {
                **series[key],
                "recipe_core_ids": sorted(series[key]["recipe_core_ids"]),
                "canonical_reaction_ids": sorted(
                    series[key]["canonical_reaction_ids"]
                ),
            }
            for key in sorted(series)
        ),
    )
    return {
        name: {
            "path": str(path),
            "sha256": _sha256(destination / path),
            "row_count": sum(1 for _ in iter_gzip_jsonl(destination / path)),
        }
        for name, path in paths.items()
    }


def _merge_shards(manifest: Mapping[str, Any], destination: Path) -> Dict[str, Any]:
    output = destination / "records.jsonl.gz"
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("wb") as target:
        for entry in manifest["shards"]:
            with (destination / entry["output_path"]).open("rb") as source:
                shutil.copyfileobj(source, target)
    temporary.replace(output)
    return {
        "path": output.name,
        "sha256": _sha256(output),
        "row_count": sum(int(entry["output_row_count"]) for entry in manifest["shards"]),
    }


def validate_sharded_conversion(
    manifest_path: str | Path,
    *,
    verify_rows: bool = True,
) -> Dict[str, Any]:
    """Validate shard provenance, checksums, schemas, and row identities."""
    source = Path(manifest_path)
    manifest = json.loads(source.read_text(encoding="utf-8"))
    destination = source.parent
    issues = []
    if manifest.get("schema_version") != SHARD_MANIFEST_SCHEMA_VERSION:
        issues.append("unsupported_manifest_schema")
    stored_contract = manifest.get("definition_contract") or {}
    stored_mapping = stored_contract.get("external_atom_mapping") or {}
    validation_provider = (
        RxnMapperProvider() if stored_mapping.get("enabled") else None
    )
    if stored_contract != _definition_contract(validation_provider):
        issues.append("stale_definition_contract")
    observation_ids = set()
    duplicate_observations = 0
    verified_rows = 0
    covered_source_rows = Counter()
    for entry in manifest.get("shards") or ():
        covered_source_rows[str(entry["source_path"])] += int(
            entry["input_row_count"]
        )
        if entry.get("status") != "complete":
            issues.append(f"incomplete_shard:{entry.get('shard_id')}")
            continue
        source_path = Path(entry["source_path"])
        if not source_path.is_file() or _sha256(source_path) != entry.get(
            "source_sha256"
        ):
            issues.append(f"source_checksum_mismatch:{entry['shard_id']}")
        output = destination / entry["output_path"]
        if not output.is_file() or _sha256(output) != entry.get("output_sha256"):
            issues.append(f"output_checksum_mismatch:{entry['shard_id']}")
            continue
        if verify_rows:
            row_count = 0
            for record in iter_gzip_jsonl(output):
                row_count += 1
                if record.get("schema_version") != RECOMMENDATION_RECORD_SCHEMA_VERSION:
                    issues.append(f"record_schema_mismatch:{entry['shard_id']}")
                if record.get("converter_definition_version") != (
                    GENERIC_CONVERTER_DEFINITION_VERSION
                ):
                    issues.append(f"converter_version_mismatch:{entry['shard_id']}")
                observation_id = str(record.get("observation_id") or "")
                if observation_id in observation_ids:
                    duplicate_observations += 1
                elif observation_id:
                    observation_ids.add(observation_id)
            verified_rows += row_count
            if row_count != int(entry["output_row_count"]):
                issues.append(f"row_count_mismatch:{entry['shard_id']}")
    for source_entry in manifest.get("source_files") or ():
        source_path = str(source_entry["path"])
        if covered_source_rows[source_path] != int(
            source_entry["covered_row_count"]
        ):
            issues.append(f"source_coverage_mismatch:{source_path}")
        if source_entry.get("coverage_complete"):
            actual_rows = sum(1 for _ in iter_csv_records(source_path))
            if actual_rows != int(source_entry["covered_row_count"]):
                issues.append(f"incomplete_source_coverage:{source_path}")
    report = {
        "schema_version": "1.0",
        "artifact_type": "sharded_conversion_integrity",
        "manifest_path": str(source),
        "valid": not issues and duplicate_observations == 0,
        "issues": sorted(set(issues)),
        "shard_count": len(manifest.get("shards") or ()),
        "verified_row_count": verified_rows if verify_rows else None,
        "duplicate_observation_count": duplicate_observations,
    }
    _atomic_json(destination / "integrity_report.json", report)
    return report


def convert_datasets_sharded(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    shard_size: int = 1_000,
    max_shards: Optional[int] = None,
    mode: str = "full",
    workers: int = 1,
    checkpoint_interval: int = 10,
    merge_records: bool = True,
    use_rxnmapper: bool = False,
    progress: bool = False,
    progress_callback: Optional[
        Callable[[ShardedConversionProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Convert a corpus into deterministic, restartable canonical shards."""
    if shard_size < 1:
        raise ValueError("shard_size must be positive")
    if max_shards is not None and max_shards < 1:
        raise ValueError("max_shards must be positive")
    if workers < 1:
        raise ValueError("workers must be positive")
    if use_rxnmapper and workers != 1:
        raise ValueError(
            "RXNMapper sharded conversion requires workers=1 to avoid loading "
            "one mapper model per process"
        )
    if checkpoint_interval < 1:
        raise ValueError("checkpoint_interval must be positive")
    paths = discover_csv_datasets(dataset_path)
    if not paths:
        raise ValueError(f"No CSV datasets found at {dataset_path}")

    def notify(
        phase: str,
        message: str,
        *,
        shard_count: int = 0,
        row_count: int = 0,
    ) -> None:
        if progress_callback is not None:
            progress_callback(
                ShardedConversionProgress(
                    phase=phase,
                    source_file_count=len(paths),
                    shard_count=shard_count,
                    row_count=row_count,
                    message=message,
                )
            )

    notify(
        "discovered",
        f"Found {len(paths)} CSV dataset file(s).",
    )
    destination = Path(output_dir)
    shards_dir = destination / "shards"
    shards_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = destination / "shard_manifest.json"
    previous_entries = {}
    if manifest_path.is_file():
        previous_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous_manifest.get("schema_version") != SHARD_MANIFEST_SCHEMA_VERSION:
            raise ValueError("Unsupported prior shard manifest schema")
        previous_entries = {
            str(entry["shard_id"]): entry
            for entry in previous_manifest.get("shards") or ()
        }
    mapping_provider = None
    if use_rxnmapper:
        if not RxnMapperProvider.is_available():
            raise RuntimeError(
                "RXNMapper is not installed; run "
                "'python -m pip install -r requirements-mapping.txt'"
            )
        mapping_provider = RxnMapperProvider()
    contract = _definition_contract(mapping_provider)
    entries = []
    source_row_counts = Counter()
    source_checksums = {}
    for file_number, path in enumerate(paths, start=1):
        if cancel_check is not None and cancel_check():
            raise ShardedConversionCancelled(
                "Conversion cancelled before shard processing"
            )
        notify(
            "hashing",
            f"Checking source file {file_number}/{len(paths)}: {path.name}",
        )
        source_checksums[str(path.resolve())] = _sha256(path)
    processed_shards = 0
    accepted_since_checkpoint = 0
    accepted_row_count = 0
    stop = False
    cancelled = False

    def checkpoint(*, coverage_complete: bool = False) -> None:
        _atomic_json(
            manifest_path,
            _manifest_payload(
                dataset_path=dataset_path,
                mode=mode,
                shard_size=shard_size,
                max_shards=max_shards,
                contract=contract,
                paths=paths,
                source_checksums=source_checksums,
                source_row_counts=source_row_counts,
                entries=entries,
                coverage_complete=coverage_complete,
            ),
        )

    def accept(entry: Mapping[str, Any]) -> None:
        nonlocal accepted_row_count, accepted_since_checkpoint
        entries.append(dict(entry))
        accepted_since_checkpoint += 1
        accepted_row_count += int(entry.get("input_row_count") or 0)
        if accepted_since_checkpoint >= checkpoint_interval:
            checkpoint()
            accepted_since_checkpoint = 0
        if progress:
            print(
                json.dumps(
                    {
                        "shard_id": entry["shard_id"],
                        "status": entry["status"],
                        "reused": entry["reused"],
                        "rows": entry["input_row_count"],
                    }
                ),
                flush=True,
            )
        notify(
            "shard_completed",
            (
                f"{'Reused' if entry.get('reused') else 'Converted' if entry.get('status') == 'complete' else 'Failed'} "
                f"{entry['shard_id']}: {entry['input_row_count']} row(s)."
            ),
            shard_count=len(entries),
            row_count=accepted_row_count,
        )

    executor = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
    shared_cache = GenericConversionCache() if executor is None else None
    pending = set()
    try:
        for file_number, path in enumerate(paths, start=1):
            source_path = str(path.resolve())
            source_sha256 = source_checksums[source_path]
            notify(
                "source_processing",
                (
                    f"Processing {file_number}/{len(paths)} (total): "
                    f"{path.name}"
                ),
                shard_count=len(entries),
                row_count=accepted_row_count,
            )
            for part_number, raw_records in enumerate(
                _chunks(iter_csv_records(path), shard_size)
            ):
                if cancel_check is not None and cancel_check():
                    cancelled = True
                    stop = True
                    break
                if max_shards is not None and processed_shards >= max_shards:
                    stop = True
                    break
                shard_id = _shard_id(path, part_number)
                output_path = Path("shards") / f"{shard_id}.jsonl.gz"
                expected = {
                    "shard_id": shard_id,
                    "source_path": source_path,
                    "source_sha256": source_sha256,
                    "part_number": part_number,
                    "input_position_start": part_number * shard_size,
                    "input_position_end": (
                        part_number * shard_size + len(raw_records) - 1
                    ),
                    "source_row_number_start": raw_records[0].source_row_number,
                    "source_row_number_end": raw_records[-1].source_row_number,
                    "input_row_count": len(raw_records),
                    "definition_contract": contract,
                    "output_path": str(output_path),
                    "mode": mode,
                }
                source_row_counts[source_path] += len(raw_records)
                processed_shards += 1
                previous = previous_entries.get(shard_id) or {}
                if _reusable(previous, expected, destination):
                    accept({**previous, "reused": True})
                    continue
                task = (
                    expected,
                    raw_records,
                    str((destination / output_path).resolve()),
                )
                if executor is None:
                    accept(
                        _convert_shard_task(
                            task,
                            cache=shared_cache,
                            mapping_provider=mapping_provider,
                        )
                    )
                    continue
                pending.add(executor.submit(_convert_shard_task, task))
                if len(pending) >= workers * 2:
                    done, pending = wait(
                        pending,
                        return_when=FIRST_COMPLETED,
                    )
                    for future in done:
                        accept(future.result())
            if stop:
                break
        if pending:
            done, _ = wait(pending)
            for future in done:
                accept(future.result())
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=False)
    if cancelled:
        checkpoint()
        notify(
            "cancelled",
            (
                f"Cancellation saved {len(entries)} completed shard(s); "
                "run again with the same output folder to resume."
            ),
            shard_count=len(entries),
            row_count=accepted_row_count,
        )
        raise ShardedConversionCancelled(
            "Conversion cancelled; completed shards were saved for resume"
        )
    manifest = _manifest_payload(
        dataset_path=dataset_path,
        mode=mode,
        shard_size=shard_size,
        max_shards=max_shards,
        contract=contract,
        paths=paths,
        source_checksums=source_checksums,
        source_row_counts=source_row_counts,
        entries=entries,
        coverage_complete=max_shards is None,
    )
    _atomic_json(manifest_path, manifest)
    complete_entries = [
        entry for entry in manifest["shards"] if entry["status"] == "complete"
    ]
    failed_entries = [
        entry for entry in manifest["shards"] if entry["status"] != "complete"
    ]
    catalogs = {}
    merged = None
    if not failed_entries:
        notify(
            "catalogs",
            "Building compressed recipe and reference catalogs…",
            shard_count=len(entries),
            row_count=sum(
                int(entry.get("input_row_count") or 0) for entry in entries
            ),
        )
        catalogs = _write_catalogs(manifest, destination)
        if merge_records:
            notify(
                "merging",
                "Merging canonical shards into records.jsonl.gz…",
                shard_count=len(entries),
                row_count=sum(
                    int(entry.get("input_row_count") or 0) for entry in entries
                ),
            )
            merged = _merge_shards(manifest, destination)
        else:
            previous_merged = destination / "records.jsonl.gz"
            if previous_merged.is_file():
                previous_merged.unlink()
    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "artifact_type": "generic_sharded_conversion_report",
        "manifest_path": str(manifest_path),
        "mode": mode,
        "merge_records": merge_records,
        "source_file_count": len(paths),
        "shard_count": len(entries),
        "complete_shard_count": len(complete_entries),
        "failed_shard_count": len(failed_entries),
        "failed_shards": [
            {
                "shard_id": str(entry.get("shard_id") or ""),
                "source_path": str(entry.get("source_path") or ""),
                "failures": list(entry.get("failures") or ()),
            }
            for entry in failed_entries
        ],
        "reused_shard_count": sum(int(entry.get("reused", False)) for entry in entries),
        "input_row_count": sum(int(entry["input_row_count"]) for entry in entries),
        "output_row_count": sum(
            int(entry.get("output_row_count") or 0) for entry in entries
        ),
        "admission_tier_counts": _merge_counts(
            complete_entries, "admission_tier_counts"
        ),
        "chemistry_status_counts": _merge_counts(
            complete_entries, "chemistry_status_counts"
        ),
        "condition_status_counts": _merge_counts(
            complete_entries, "condition_status_counts"
        ),
        "condition_stage_status_counts": _merge_counts(
            complete_entries, "condition_stage_status_counts"
        ),
        "outcome_status_counts": _merge_counts(
            complete_entries, "outcome_status_counts"
        ),
        "index_eligibility_counts": _merge_counts(
            complete_entries, "index_eligibility_counts"
        ),
        "admission_reason_counts": _merge_counts(
            complete_entries, "admission_reason_counts"
        ),
        "transformation_class_counts": _merge_counts(
            complete_entries, "transformation_class_counts"
        ),
        "named_family_counts": _merge_counts(
            complete_entries, "named_family_counts"
        ),
        "reaction_scope_counts": _merge_counts(
            complete_entries, "reaction_scope_counts"
        ),
        "reaction_completeness_status_counts": _merge_counts(
            complete_entries, "reaction_completeness_status_counts"
        ),
        "fragment_source_support_status_counts": _merge_counts(
            complete_entries,
            "fragment_source_support_status_counts",
        ),
        "external_atom_mapping": {
            **contract["external_atom_mapping"],
            "status_counts": _merge_counts(
                complete_entries,
                "external_mapping_status_counts",
            ),
        },
        "signature_count": sum(
            int(entry.get("signature_count") or 0) for entry in complete_entries
        ),
        "catalogs": catalogs,
        "merged_records": merged,
    }
    _atomic_json(destination / "conversion_report.json", report)
    integrity = validate_sharded_conversion(
        manifest_path,
        verify_rows=not failed_entries,
    )
    report["integrity"] = integrity
    _atomic_json(destination / "conversion_report.json", report)
    notify(
        "completed",
        (
            f"Canonical conversion complete: "
            f"{report['output_row_count']} record(s)."
        ),
        shard_count=len(entries),
        row_count=int(report["output_row_count"]),
    )
    return report


__all__ = [
    "SHARDED_CONVERSION_DEFINITION_VERSION",
    "SHARD_MANIFEST_SCHEMA_VERSION",
    "ShardedConversionCancelled",
    "ShardedConversionProgress",
    "convert_datasets_sharded",
    "iter_gzip_jsonl",
    "validate_sharded_conversion",
]
