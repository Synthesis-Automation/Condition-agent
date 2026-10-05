"""Resumable preparation of a single complete processed dataset release."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

from .conversion.atomic import atomic_json
from .conversion.input_schema import discover_conversion_datasets
from .conversion.sharded import _definition_contract, convert_datasets_sharded
from .corpus_io import file_sha256
from .record_storage import STORAGE_SCHEMA_VERSION

DEFAULT_INPUT_ROOT = Path("datasets/intermediate_datasets")
DEFAULT_OUTPUT_ROOT = Path("datasets/processed_datasets")
RELEASE_SCHEMA_VERSION = "processed_reaction_release.v1"


def prepare_release(
    source: str | Path, output_root: str | Path,
) -> tuple[Path, dict[str, Any]]:
    """Fingerprint all declared observations and choose a reproducible build root."""
    root = Path(source).resolve()
    sources = discover_conversion_datasets(root)
    if not sources:
        raise ValueError("No physical reaction observations found")
    inventory = [{"path": str(p.resolve()), "sha256": file_sha256(p),
                  "size_bytes": p.stat().st_size} for p in sources]
    route_sources = sorted(root.rglob("routes.source.jsonl.gz")) if root.is_dir() else []
    route_inventory = [{"path": str(p.resolve()), "sha256": file_sha256(p),
                        "size_bytes": p.stat().st_size} for p in route_sources]
    contract = _definition_contract()
    repository = Path(__file__).resolve().parents[1]
    definition_files = {path.relative_to(repository).as_posix(): file_sha256(path)
                        for package in ("reactive_taxonomy", "condition_registry", "condition_recommender",
                                        "core_retrosynthesis", "forward_synthesis")
                        for path in sorted((repository / package).rglob("*.json"))
                        if "definitions" in path.parts}
    identity = {"sources": [entry["sha256"] for entry in inventory],
                "definition_contract": contract, "storage_schema": STORAGE_SCHEMA_VERSION,
                "definition_files": definition_files,
                "route_sources": [entry["sha256"] for entry in route_inventory],
                "pipeline_generation": "processed_pipeline.v3", "release_schema": RELEASE_SCHEMA_VERSION}
    release_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:24]
    release = Path(output_root).resolve() / "releases" / release_id
    release.mkdir(parents=True, exist_ok=True)
    manifest = {"schema_version": RELEASE_SCHEMA_VERSION, "release_id": release_id,
                "build_complete": False, "source_coverage": "all_physical_observations",
                "sources": inventory, "definition_contract": contract,
                "definition_files": definition_files,
                "route_sources": route_inventory, "pipeline_generation": "processed_pipeline.v3",
                "storage_schema": STORAGE_SCHEMA_VERSION,
                "independent_chemistry_review": "pending",
                "conversion_authorization": "user_requested_full_corpus_conversion_2026-10-05"}
    return release, manifest


def convert_release(
    source: str | Path = DEFAULT_INPUT_ROOT,
    output_root: str | Path = DEFAULT_OUTPUT_ROOT, *, workers: int = 1,
    shard_size: int = 500, progress: Callable[[dict[str, Any]], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Convert every observation into lossless shared-object canonical shards."""
    release, manifest = prepare_release(source, output_root)
    if not (release / "manifest.json").exists():
        atomic_json(release / "manifest.json", manifest)
    def emit(event: Any) -> None:
        if progress:
            progress({"phase": event.phase, "rows": event.row_count,
                      "shards": event.shard_count, "message": event.message})
    report = convert_datasets_sharded(source, release / "records", workers=workers,
                                     shard_size=shard_size, merge_records=False,
                                     build_catalogs=False, progress_callback=emit,
                                     cancel_check=cancel_check)
    if report["failed_shard_count"] or not report["integrity"]["valid"]:
        raise RuntimeError("Canonical conversion is incomplete; inspect conversion_report.json")
    return {"release_path": str(release), "release_id": manifest["release_id"],
            "conversion": report}


def main() -> None:
    """Build all observations without a production sampling mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("convert",))
    parser.add_argument("--source", default=str(DEFAULT_INPUT_ROOT))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shard-size", type=int, default=500)
    args = parser.parse_args()
    result = convert_release(args.source, args.output_root, workers=args.workers,
                             shard_size=args.shard_size,
                             progress=lambda event: print(json.dumps(event), flush=True))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
