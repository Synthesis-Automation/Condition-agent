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
from .processed_catalog import CATALOG_SCHEMA_VERSION
from .record_storage import STORAGE_SCHEMA_VERSION

DEFAULT_INPUT_ROOT = Path("datasets/intermediate_datasets")
DEFAULT_OUTPUT_ROOT = Path("datasets/processed_datasets")
RELEASE_SCHEMA_VERSION = "processed_reaction_release.v1"
PIPELINE_GENERATION = "processed_pipeline.v5"
ROUTE_ARCHIVE = "routes.source.jsonl.gz"
ROUTE_STEPS = "route_steps.observations.jsonl.gz"


def discover_route_bundles(source: str | Path) -> list[tuple[Path, Path]]:
    """Resolve declared route/observation pairs, including older migration inventories."""
    root = Path(source).resolve()
    if not root.is_dir():
        return []
    archives = set(root.rglob(ROUTE_ARCHIVE))
    for release_manifest in root.rglob("route_release.manifest.json"):
        if release_manifest.with_name(ROUTE_ARCHIVE) not in archives:
            raise ValueError("Missing route archive declared by route release manifest")
    inventory = root / "migration.manifest.json"
    if inventory.exists():
        manifest = json.loads(inventory.read_text(encoding="utf-8"))
        if not manifest.get("coverage_complete"):
            raise ValueError("Migration source coverage is incomplete")
        declared = []
        total = 0
        for entry in manifest["sources"]:
            if entry["category"] != "routes_and_single_step":
                continue
            relative = entry.get("bundle_relative_path") or Path(entry["source_relative_path"]).parent.parent
            folder = (root / relative).resolve()
            if not folder.is_relative_to(root):
                raise ValueError("Route bundle escapes its intermediate root")
            declared.append(folder / ROUTE_ARCHIVE)
            total += entry["report"]["counts"]["route_count"]
        if (len(declared) != len(set(declared)) or set(declared) != archives
            or total != manifest["route_count"] or (total > 0 and not archives)):
            raise ValueError("Declared route bundles do not match migration inventory")
    pairs = [(archive, archive.with_name(ROUTE_STEPS)) for archive in sorted(archives)]
    if any(not steps.is_file() for _, steps in pairs):
        raise ValueError("Missing observation archive for a route bundle")
    return pairs


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
    pairs = discover_route_bundles(root)
    observation_hashes = {row["path"]: row["sha256"] for row in inventory}
    route_inventory = []
    for routes, steps in pairs:
        item = {"path": str(routes), "sha256": file_sha256(routes),
                "size_bytes": routes.stat().st_size, "step_path": str(steps),
                "step_sha256": observation_hashes[str(steps)]}
        release_manifest = routes.with_name("route_release.manifest.json")
        if release_manifest.is_file():
            declared = json.loads(release_manifest.read_text(encoding="utf-8"))
            hashes = declared["output_sha256"]
            if (not declared["coverage_complete"] or hashes[ROUTE_ARCHIVE] != item["sha256"]
                or hashes[ROUTE_STEPS] != item["step_sha256"]):
                raise ValueError(f"Route bundle differs from its complete manifest: {routes.parent}")
            item.update(expected_route_count=declared["counts"]["route_count"],
                        expected_observation_count=declared["counts"]["output_step_observations"])
        route_inventory.append(item)
    contract = _definition_contract()
    repository = Path(__file__).resolve().parents[1]
    definition_files = {path.relative_to(repository).as_posix(): file_sha256(path)
                        for package in ("reactive_taxonomy", "condition_registry", "condition_recommender",
                                        "core_retrosynthesis", "forward_synthesis")
                        for path in sorted((repository / package).rglob("*.json"))
                        if "definitions" in path.parts}
    identity = {"sources": [entry["sha256"] for entry in inventory],
                "definition_contract": contract, "storage_schema": STORAGE_SCHEMA_VERSION,
                "evidence_catalog_schema": CATALOG_SCHEMA_VERSION,
                "definition_files": definition_files,
                "route_sources": [(e["sha256"], e["step_sha256"], e.get("expected_route_count"),
                                   e.get("expected_observation_count")) for e in route_inventory],
                "pipeline_generation": PIPELINE_GENERATION, "release_schema": RELEASE_SCHEMA_VERSION}
    release_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:24]
    release = Path(output_root).resolve() / "releases" / release_id
    release.mkdir(parents=True, exist_ok=True)
    manifest = {"schema_version": RELEASE_SCHEMA_VERSION, "release_id": release_id,
                "build_complete": False, "source_coverage": "all_physical_observations",
                "sources": inventory, "definition_contract": contract,
                "definition_files": definition_files,
                "route_sources": route_inventory, "pipeline_generation": PIPELINE_GENERATION,
                "storage_schema": STORAGE_SCHEMA_VERSION,
                "evidence_catalog_schema": CATALOG_SCHEMA_VERSION,
                "independent_chemistry_review": "pending",
                "conversion_authorization": "requested_build"}
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
