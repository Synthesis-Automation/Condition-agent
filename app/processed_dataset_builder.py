"""Build and atomically publish all app/agent artifacts from one complete corpus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import monotonic
from typing import Any, Callable

from condition_recommender.conversion.atomic import atomic_json
from condition_recommender.conversion.sharded import validate_sharded_conversion
from condition_recommender.corpus_io import file_sha256, iter_canonical_records
from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.generic_indexing import load_generic_index
from condition_recommender.processed_build import (
    DEFAULT_INPUT_ROOT, DEFAULT_OUTPUT_ROOT, PIPELINE_GENERATION, convert_release,
    observation_source_binding, prepare_release, verify_route_bundle_report,
)
from condition_recommender.processed_catalog import CATALOG_SCHEMA_VERSION, build_processed_catalog
from condition_recommender.processed_release import publish_processed_release, resolve_processed_release
from condition_recommender.shared_core_index import build_shared_core_index, load_shared_core_index
from condition_recommender.sqlite_indexing import build_sqlite_generic_index
from core_retrosynthesis.full_scale import build_full_scale_operator_library
from core_retrosynthesis.processed_routes import ROUTE_CATALOG_VERSION, build_route_bundle_catalog
from core_retrosynthesis.route_curation import ROUTE_CURATION_POLICY_VERSION
from forward_synthesis import build_forward_library, save_forward_library
from reactive_taxonomy.shared_core_graph import shared_core_availability


def build_processed_datasets(
    source: str | Path = DEFAULT_INPUT_ROOT,
    output_root: str | Path = DEFAULT_OUTPUT_ROOT, *, workers: int = 1,
    shard_size: int = 500, finish_only: bool = False,
    progress: Callable[[dict[str, Any]], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Build all required artifacts, then publish their coherent release manifest."""
    release, manifest = prepare_release(source, output_root)
    reports = release / "reports"
    reports.mkdir(exist_ok=True)
    state_path = reports / "build_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {}
    if (release / "manifest.json").exists():
        previous = json.loads((release / "manifest.json").read_text())
        if previous.get("build_complete") and previous.get("builder_version") == PIPELINE_GENERATION:
            return publish_processed_release(release, output_root).manifest

    def emit(phase: str, **values: Any) -> None:
        if cancel_check and cancel_check():
            raise InterruptedError("Build cancelled; completed stages remain resumable")
        if progress:
            progress({"phase": phase, **values})

    if not finish_only:
        convert_release(source, output_root, workers=workers, shard_size=shard_size,
                        progress=progress, cancel_check=cancel_check)
    canonical = release / "records" / "shard_manifest.json"
    if not canonical.exists():
        raise ValueError("Canonical conversion has not run")
    conversion = json.loads(canonical.read_text())
    if (conversion.get("coverage_complete") is False
        or not conversion.get("source_files")
        or any(not s.get("coverage_complete") for s in conversion["source_files"])
        or any(s.get("status") != "complete" for s in conversion["shards"])):
        raise ValueError("Canonical conversion is unfinished")
    if observation_source_binding(conversion["source_files"]) != observation_source_binding(manifest["sources"]):
        raise ValueError("Canonical source hashes differ from the selected release inventory")
    binding = file_sha256(canonical)
    artifacts: dict[str, dict[str, Any]] = {}

    def entry(path: Path) -> dict[str, Any]:
        return {"relative_path": path.relative_to(release).as_posix(),
                "size_bytes": path.stat().st_size, "sha256": file_sha256(path)}

    def stage(name: str, outputs: dict[str, Path], run: Callable[[], Any], version: str,
              dependencies: tuple[str, ...] = ()) -> Any:
        emit(name, message=f"Building {name}")
        started = monotonic()
        old = state.get(name) or {}
        dependency_bindings = {key: artifacts[key]["sha256"] for key in dependencies}
        if old.get("canonical_binding") == binding and old.get("builder_version") == version and old.get("dependency_bindings", {}) == dependency_bindings and all(
            p.is_file() and file_sha256(p) == old.get("artifacts", {}).get(k, {}).get("sha256")
            for k, p in outputs.items()
        ):
            result = old.get("report")
            emit(name + "_reused")
        else:
            result = run()
        current = {k: entry(p) for k, p in outputs.items()}
        artifacts.update(current)
        state[name] = {"canonical_binding": binding, "builder_version": version,
                       "dependency_bindings": dependency_bindings,
                       "artifacts": current, "report": result,
                       "elapsed_seconds": round(monotonic() - started, 3)}
        atomic_json(state_path, state)
        emit(name + "_completed", report=result)
        return result

    indexes = release / "indexes"
    indexes.mkdir(exist_ok=True)
    availability = shared_core_availability().fingerprint
    catalog_report = stage("catalogs", {"evidence_catalog": indexes / "catalogs.sqlite"},
          lambda: build_processed_catalog(canonical, indexes / "catalogs.sqlite",
                      progress=lambda event: emit("catalog_rows", report=event)), f"{CATALOG_SCHEMA_VERSION}/{availability}")
    route_sources = manifest["route_sources"]
    route_report = None
    if route_sources:
        pairs = [(item["path"], item["step_path"]) for item in route_sources]
        route_report = stage("routes", {"route_catalog": indexes / "route_catalog.sqlite"},
              lambda: build_route_bundle_catalog(pairs, indexes / "route_catalog.sqlite",
                      workers=workers, evidence_catalog=indexes / "catalogs.sqlite",
                      progress=lambda event: emit("route_rows", report=event)),
                      f"{ROUTE_CATALOG_VERSION}/curation.v{ROUTE_CURATION_POLICY_VERSION}",
                      ("evidence_catalog",))
        for source_entry, result in zip(route_sources, route_report["bundles"], strict=True):
            verify_route_bundle_report(source_entry, result)
    condition_report = stage("conditions", {"condition_index": indexes / "generic_index.sqlite"},
          lambda: build_sqlite_generic_index(iter_canonical_records(canonical), indexes / "generic_index.sqlite",
                                              progress_callback=lambda phase, count: emit("condition_" + phase, rows=count),
                                              cancel_check=cancel_check), "generic_index.v6.5")
    index = load_generic_index(indexes / "generic_index.sqlite")
    stage("shared_core", {"shared_core_index": indexes / "generic_index.shared_core.sqlite"},
          lambda: build_shared_core_index(index, indexes / "generic_index.shared_core.sqlite",
                                          workers=workers, resume=True,
                                          progress_callback=lambda count: emit("shared_core_rows", rows=count),
                                          cancel_check=cancel_check), f"shared_core.v3/{availability}", ("condition_index",))
    fragment_report = stage("fragments", {"fragment_index": indexes / "fragment_index.sqlite"},
          lambda: build_fragment_index(canonical, indexes / "fragment_index.sqlite",
                                        evidence_catalog=indexes / "catalogs.sqlite",
                                        progress=lambda event: emit("fragment_rows", **event)), "fragments.v4", ("evidence_catalog",))

    def operators() -> dict[str, Any]:
        library, report = build_full_scale_operator_library(
            canonical, release / "operators" / "retrosynthesis", workers=workers,
            progress_callback=lambda event: emit("retro_" + event.get("phase", "build"), report=event))
        forward = build_forward_library(library, workers=workers)
        save_forward_library(forward, release / "operators" / "forward" / "forward_operator_library_v1.json.gz")
        return report

    stage("operators", {
        "retro_library": release / "operators" / "retrosynthesis" / "operator_library_v3.json.gz",
        "forward_library": release / "operators" / "forward" / "forward_operator_library_v1.json.gz",
    }, operators, "operators.v4")

    load_shared_core_index(indexes / "generic_index.shared_core.sqlite", index)
    emit("validation", message="Checking complete source coverage and all canonical shard identities")
    integrity = validate_sharded_conversion(canonical)
    if not integrity["valid"]:
        raise ValueError("Canonical validation failed")
    total = sum(s["output_row_count"] for s in conversion["shards"])
    if catalog_report["counts"]["observations"] != total or fragment_report["counts"]["source_observations"] != total:
        raise ValueError("Derived evidence/fragment coverage does not match the complete corpus")
    expected_path = Path(source) / "migration.manifest.json"
    if expected_path.exists():
        expected_manifest = json.loads(expected_path.read_text())
        expected = expected_manifest["single_step_observation_count"]
        if total != expected:
            raise ValueError(f"Expected {expected} physical observations; converted {total}")
        if (route_report or {}).get("route_count", 0) != expected_manifest["route_count"]:
            raise ValueError("Source route count does not match migration inventory")
    coverage = {"observation_count": total, "source_count": len(conversion["source_files"]),
                "source_coverage_complete": True, "validation": integrity,
                "routes": route_report,
                "evidence_catalog": catalog_report, "condition_index": condition_report,
                "fragment_index": fragment_report,
                "legacy_organic_syntheses_example": "retired_from_production_20_row_example_not_in_declared_raw_corpus",
                "algorithmic_abstractions": "excluded_not_physical_reaction_observations",
                "scientific_review_status": "independent_review_pending"}
    atomic_json(reports / "coverage.json", coverage)
    artifacts["canonical_records"] = entry(canonical)
    artifacts["coverage_report"] = entry(reports / "coverage.json")
    manifest.update(build_complete=True, builder_version=PIPELINE_GENERATION, artifacts=artifacts, coverage=coverage,
                    capabilities={**{name: True for name in artifacts}, "composite_promotion": False},
                    capability_limits={"composite_promotion": "Requires declared training scope and independent support review; source routes remain accessible."})
    atomic_json(release / "manifest.json", manifest)
    publish_processed_release(release, output_root)
    emit("published", release_id=manifest["release_id"], observations=total)
    return manifest


def main() -> None:
    """Run a complete resumable build or finish previously converted canonical shards."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "finish", "validate"))
    parser.add_argument("--source", default=str(DEFAULT_INPUT_ROOT))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shard-size", type=int, default=500)
    args = parser.parse_args()
    if args.command == "validate":
        result = resolve_processed_release(args.output_root, verify_artifacts=True).manifest
    else:
        result = build_processed_datasets(args.source, args.output_root, workers=args.workers,
                                          shard_size=args.shard_size, finish_only=args.command == "finish",
                                          progress=lambda event: print(json.dumps(event), flush=True))
    print(json.dumps({"release_id": result["release_id"], "coverage": result["coverage"]}), flush=True)


if __name__ == "__main__":
    main()
