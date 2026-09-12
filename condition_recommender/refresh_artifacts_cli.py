"""Restartable artifact refresh from source-faithful canonical observation shards."""

from __future__ import annotations

import argparse
import json
import shutil
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from rdkit import RDLogger

from .conversion.atomic import atomic_json, atomic_output_path
from .conversion.generic import GenericConversionCache
from .conversion.refresh import refresh_canonical_record
from .conversion.sharded import (
    _converted_counts,
    _definition_contract,
    _sha256,
    _write_gzip_jsonl,
    iter_gzip_jsonl,
    validate_sharded_conversion,
    write_conversion_catalogs,
)
from .sqlite_indexing import build_sqlite_generic_index


def _refresh_shard(task: dict[str, Any]) -> dict[str, Any]:
    RDLogger.DisableLog("rdApp.*")
    output = Path(task["output"])
    output.parent.mkdir(parents=True, exist_ok=True)
    checkpoint = output.with_suffix(".refresh.json")
    identity = {key: task[key] for key in ("input", "input_sha256", "contract")}
    identity["refresh_algorithm_version"] = "1.0"
    if checkpoint.is_file() and output.is_file():
        saved = json.loads(checkpoint.read_text(encoding="utf-8"))
        if (
            saved.get("identity") == identity
            and _sha256(output) == saved["entry"]["output_sha256"]
        ):
            return saved["entry"]
    source = Path(task["input"])
    if _sha256(source) != task["input_sha256"]:
        raise ValueError(f"Source shard checksum mismatch: {source}")
    cache = GenericConversionCache(max_entries=256)
    records = []
    changes = Counter()
    started = time.perf_counter()
    for old in iter_gzip_jsonl(source):
        new = refresh_canonical_record(old, cache=cache)
        for field in (
            "admission_tier",
            "chemistry_status",
            "condition_status",
            "index_eligibility",
            "precedent_tier",
            "named_family",
        ):
            if old.get(field) != new.get(field):
                changes[f"{field}:{old.get(field)}->{new.get(field)}"] += 1
        records.append(new)
    count = _write_gzip_jsonl(output, records)
    if count != int(task["entry"]["output_row_count"]):
        raise ValueError("Canonical refresh changed source row count")
    entry = {
        **task["entry"],
        **_converted_counts(records),
        "definition_contract": task["contract"],
        "output_path": task["output_path"],
        "output_sha256": _sha256(output),
        "output_row_count": count,
        "status": "complete",
        "failure_count": 0,
        "failures": [],
        "reused": False,
        "refresh_source_path": str(source),
        "refresh_source_sha256": task["input_sha256"],
        "refresh_changes": dict(changes),
        "refresh_seconds": round(time.perf_counter() - started, 3),
    }
    atomic_json(checkpoint, {"identity": identity, "entry": entry})
    return entry


def refresh_saved_library(
    source_library: Path,
    output: Path,
    *,
    workers: int = 8,
    max_shards: int | None = None,
    build_index: bool = True,
) -> dict[str, Any]:
    """Refresh an inventoried saved batch into a separate, integrity-checked library."""
    source_library, output = source_library.resolve(), output.resolve()
    if source_library == output or source_library in output.parents:
        raise ValueError(
            "Refresh output must be separate from the frozen source library"
        )
    if workers < 1:
        raise ValueError("workers must be positive")
    combined_path = source_library / "combined_recommendation_report.json"
    combined = json.loads(combined_path.read_text(encoding="utf-8"))
    contract = _definition_contract()
    tasks, sources, seen = [], {}, set()
    for batch_path in combined["batch_manifest_paths"]:
        batch = json.loads(Path(batch_path).read_text(encoding="utf-8"))
        for source in batch["source_files"]:
            artifact = (
                source_library / "converted_sources" / source["source_artifact_id"]
            )
            manifest = json.loads(
                (artifact / "shard_manifest.json").read_text(encoding="utf-8")
            )
            if not all(
                item.get("coverage_complete")
                for item in manifest.get("source_files") or ()
            ) or any(item.get("status") != "complete" for item in manifest["shards"]):
                raise ValueError(
                    "Refresh requires complete source conversion manifests"
                )
            if (
                manifest["definition_contract"]["taxonomy_definition_versions"]
                != contract["taxonomy_definition_versions"]
            ):
                raise ValueError(
                    "Source observation definitions are stale; use full source conversion"
                )
            for entry in manifest["shards"]:
                input_path = (artifact / entry["output_path"]).resolve()
                if input_path in seen:
                    continue
                seen.add(input_path)
                relative = str(Path("shards") / input_path.name)
                tasks.append(
                    {
                        "input": str(input_path),
                        "input_sha256": entry["output_sha256"],
                        "output": str(output / relative),
                        "output_path": relative,
                        "entry": entry,
                        "contract": contract,
                    }
                )
            for value in manifest["source_files"]:
                sources[value["path"]] = value
    tasks.sort(
        key=lambda task: (task["entry"]["source_path"], task["entry"]["part_number"])
    )
    if len({task["output_path"] for task in tasks}) != len(tasks):
        raise ValueError("Refresh shard filenames collide")
    if sum(int(task["entry"]["output_row_count"]) for task in tasks) != int(
        combined["record_count"]
    ):
        raise ValueError(
            "Source shard inventory differs from the combined canonical corpus"
        )
    if max_shards is not None:
        tasks = tasks[:max_shards]
    output.mkdir(parents=True, exist_ok=True)
    atomic_json(
        output / "refresh_source_inventory.json",
        {
            "source_library": str(source_library),
            "source_report_sha256": _sha256(combined_path),
            "contract": contract,
            "shards": [
                {key: task[key] for key in ("input", "input_sha256", "output_path")}
                for task in tasks
            ],
        },
    )
    results = {}
    print(f"Refreshing {len(tasks)} shards with {workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        pending = {executor.submit(_refresh_shard, task): task for task in tasks}
        for future in as_completed(pending):
            entry = future.result()
            results[entry["output_path"]] = entry
            if len(results) % 5 == 0 or len(results) == len(tasks):
                print(
                    f"Refreshed {len(results)}/{len(tasks)} shards; {sum(item['output_row_count'] for item in results.values())} rows",
                    flush=True,
                )
    entries = [results[task["output_path"]] for task in tasks]
    source_counts = Counter()
    for entry in entries:
        source_counts[entry["source_path"]] += entry["input_row_count"]
    source_entries = [
        {
            **sources[path],
            "covered_row_count": count,
            "coverage_complete": max_shards is None,
        }
        for path, count in sorted(source_counts.items())
    ]
    manifest = {
        "schema_version": "1.0",
        "artifact_type": "canonical_observation_refresh",
        "definition_contract": contract,
        "coverage_complete": max_shards is None,
        "source_files": source_entries,
        "shards": entries,
    }
    manifest_path = output / "shard_manifest.json"
    atomic_json(manifest_path, manifest)
    print("Validating all refreshed source and output rows", flush=True)
    integrity = validate_sharded_conversion(manifest_path)
    if not integrity["valid"]:
        atomic_json(output / "integrity_report.json", integrity)
        raise ValueError(f"Refreshed conversion integrity failed: {integrity}")
    counts = {}
    for entry in entries:
        for key, value in entry.items():
            if key.endswith("_counts") and isinstance(value, dict):
                counts.setdefault(key, Counter()).update(value)
    report = {
        "schema_version": "1.0",
        "artifact_type": "canonical_observation_refresh",
        "definition_contract": contract,
        "output_row_count": sum(item["output_row_count"] for item in entries),
        "failed_shard_count": 0,
        "shard_count": len(entries),
        "coverage_complete": max_shards is None,
        "integrity": integrity,
        **{key: dict(value) for key, value in counts.items()},
    }
    changes = Counter()
    for entry in entries:
        changes.update(entry["refresh_changes"])
    report["changes_from_frozen_records"] = dict(changes)
    atomic_json(output / "conversion_report.json", report)

    def records():
        for entry in entries:
            yield from iter_gzip_jsonl(output / entry["output_path"])

    print("Writing recipe, reference and procedure catalogs", flush=True)
    report["catalogs"] = write_conversion_catalogs(records(), output)
    merged = output / "records.jsonl.gz"
    with atomic_output_path(merged) as temporary:
        with temporary.open("wb") as target:
            for entry in entries:
                with (output / entry["output_path"]).open("rb") as source:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
    if build_index:
        print("Building the current trusted SQLite index", flush=True)
        report["index"] = build_sqlite_generic_index(
            records(),
            output / "generic_index.sqlite",
            progress_callback=lambda phase, count: print(
                f"Index {phase}: {count}", flush=True
            )
            if count % 10000 == 0
            else None,
        )
        if report.get("precedent_tier_counts", {}).get("review_core", 0):
            report["review_index"] = build_sqlite_generic_index(
                records(), output / "generic_review_index.sqlite", include_review=True
            )
    atomic_json(output / "refresh_report.json", report)
    print(f"Refresh complete: {report['output_row_count']} rows", flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_library", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--max-shards", type=int)
    parser.add_argument("--skip-index", action="store_true")
    args = parser.parse_args()
    refresh_saved_library(
        args.source_library,
        args.output,
        workers=args.workers,
        max_shards=args.max_shards,
        build_index=not args.skip_index,
    )


if __name__ == "__main__":
    main()
