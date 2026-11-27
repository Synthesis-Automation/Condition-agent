"""Rebuild a consolidated intermediate corpus from the local raw source tree."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from ..conversion.atomic import atomic_json
from .artifacts import PreprocessingProgress, _sha256, preprocess_file
from .registry import detect_adapter
from .released_sources import released_source

MIGRATION_DEFINITION_VERSION = "intermediate_corpus_migration.v1.1"


def _inventory_counts(entries: list[dict[str, Any]]) -> dict[str, int]:
    counts = {
        "source_dataset_count": len(entries),
        "single_step_observation_count": 0,
        "route_count": 0,
        "abstraction_count": 0,
    }
    for entry in entries:
        category = entry["category"]
        if category == "routes_and_single_step":
            report_counts = entry["report"]["counts"]
            counts["route_count"] += report_counts["route_count"]
            counts["single_step_observation_count"] += report_counts[
                "output_step_observations"
            ]
        elif category in {"abstractions", "single_step"}:
            key = (
                "abstraction_count"
                if category == "abstractions"
                else "single_step_observation_count"
            )
            counts[key] += entry["report"]["output_row_count"]
        elif category != "validated_export_superseded_by_route_steps":
            raise ValueError(f"Unknown intermediate category: {category}")
    return counts


def regenerate_intermediate_datasets(
    raw_root: str | Path,
    output_root: str | Path,
    *,
    force: bool = False,
    cbyd_release: str | None = None,
    progress_callback: Callable[[PreprocessingProgress], None] | None = None,
) -> dict[str, Any]:
    """Prepare all supported sources with explicit route/export semantics.

    The original route CSV is cross-checked rather than added a second time.
    Abstracted exports use a suffix excluded by recommendation input discovery.
    Each per-file preparation can be reused after a failed or interrupted run.
    """
    source, destination = Path(raw_root).resolve(), Path(output_root).resolve()
    if not source.is_dir():
        raise ValueError(f"Raw source directory does not exist: {source}")
    if (
        destination == source
        or destination.is_relative_to(source)
        or source.is_relative_to(destination)
    ):
        raise ValueError("Raw and intermediate trees must be separate")
    csvs = sorted(source.rglob("*.csv"), key=lambda p: p.relative_to(source).as_posix())
    route_paths = sorted([*source.rglob("*.jsonl.gz"), *source.rglob("cbyd-v*.tar.gz")])
    # Preflight everything before writing a large corpus.
    adapters = {path: detect_adapter(path) for path in csvs}
    releases = []
    cbyd_versions = set()
    for path in route_paths:
        release = released_source(path, source)
        if release is None:
            raise ValueError(f"Unsupported route source: {path}")
        if release.step_export is None:
            cbyd_versions.add(path.parent.name)
            if cbyd_release and path.parent.name != cbyd_release:
                continue
        elif not release.step_export.is_file():
            raise ValueError(f"Missing original reaction export for {path}")
        releases.append(release)
    if cbyd_release and cbyd_release not in cbyd_versions:
        raise ValueError(f"CbyD release was not found: {cbyd_release}")
    if not cbyd_release and len(cbyd_versions) > 1:
        raise ValueError(
            "Multiple CbyD snapshots found; select one with --cbyd-release"
        )
    if len({r.bundle_relative_path for r in releases}) != len(releases):
        raise ValueError(
            "More than one route input targets the same intermediate bundle"
        )
    declared_bundles = {destination / r.bundle_relative_path for r in releases}
    if any(
        p.parent not in declared_bundles
        for p in destination.rglob("routes.source.jsonl.gz")
    ):
        raise ValueError(
            "Output contains a different route inventory; use a separate output root"
        )
    superseded = {r.step_export for r in releases if r.step_export}
    csv_jobs = []
    expected_steps = {
        folder / "route_steps.observations.jsonl.gz" for folder in declared_bundles
    }
    for path in csvs:
        if path in superseded:
            continue
        adapter = adapters[path]
        category = (
            "abstractions"
            if adapter.adapter_id == "higher_level_abstraction_csv.v1"
            else "single_step"
        )
        output = destination / category / path.parent.relative_to(source)
        csv_jobs.append((path, category, output))
        if category == "single_step":
            suffix = getattr(adapter, "artifact_suffix", ".observations.jsonl.gz")
            expected_steps.add(output / f"{path.name}{suffix}")
    if any(
        p not in expected_steps for p in destination.rglob("*.observations.jsonl.gz")
    ):
        raise ValueError(
            "Output contains a different observation inventory; use a separate output root"
        )
    entries = []
    for path, category, output in csv_jobs:
        report = preprocess_file(
            path,
            output,
            force=force,
            progress_callback=progress_callback,
        )
        entries.append(
            {
                "source_relative_path": path.relative_to(source).as_posix(),
                "output_relative_path": Path(report["output_path"])
                .relative_to(destination)
                .as_posix(),
                "category": category,
                "report": report,
            }
        )
    for release in releases:
        report = release.prepare(
            destination,
            force=force,
            progress_callback=progress_callback,
        )
        entries.append(
            {
                "source_relative_path": release.archive.relative_to(source).as_posix(),
                "bundle_relative_path": release.bundle_relative_path.as_posix(),
                "category": "routes_and_single_step",
                "report": report,
            }
        )
        if release.step_export:
            entries.append(
                {
                    "source_relative_path": release.step_export.relative_to(
                        source
                    ).as_posix(),
                    "category": "validated_export_superseded_by_route_steps",
                    "source_sha256": report["original_csv_sha256"],
                }
            )
    if not entries:
        raise ValueError("No supported raw source datasets found")
    report = {
        "definition_version": MIGRATION_DEFINITION_VERSION,
        "raw_root": str(source),
        "output_root": str(destination),
        **_inventory_counts(entries),
        "sources": entries,
        "coverage_complete": all(
            entry["report"]["coverage_complete"]
            for entry in entries
            if "report" in entry
        ),
        "downstream_conversion_scope": (
            "All *.observations.jsonl.gz under this output root; "
            "abstractions and raw routes are excluded."
        ),
        "migration_status": "prepared_not_chemistry_admitted",
        "legacy_intermediate_removed": False,
    }
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / "migration.manifest.json"
    if not report["coverage_complete"]:
        raise ValueError("Migration coverage is incomplete")
    atomic_json(target, report)
    return report


def verified_intermediate_artifacts(output_root: str | Path) -> dict[Path, str]:
    """Verify each declared intermediate output and return its absolute path and digest."""
    root = Path(output_root).resolve()
    report = json.loads((root / "migration.manifest.json").read_text(encoding="utf-8"))
    if not report["coverage_complete"] or any(
        report.get(k) != v for k, v in _inventory_counts(report["sources"]).items()
    ):
        raise ValueError("Migration coverage or inventory counts are inconsistent")
    verified: dict[Path, str] = {}
    expected_inputs = set()
    for entry in report["sources"]:
        item = entry.get("report")
        if item is None:
            continue
        if not item.get("coverage_complete"):
            raise ValueError("Incomplete intermediate source")
        if entry["category"] in {"single_step", "abstractions"}:
            relative = entry.get("output_relative_path")
            if relative is None:
                relative = Path(item["output_path"]).relative_to(report["output_root"])
            outputs = {root / relative: item["output_sha256"]}
        elif entry["category"] == "routes_and_single_step":
            relative = (
                entry.get("bundle_relative_path")
                or Path(entry["source_relative_path"]).parent.parent
            )
            bundle = (root / relative).resolve()
            if not bundle.is_relative_to(root / "routes"):
                raise ValueError(f"Invalid intermediate bundle: {bundle}")
            bundle_report = json.loads(
                (bundle / "route_release.manifest.json").read_text(encoding="utf-8")
            )
            hashes = item["output_sha256"]
            if set(hashes) != {
                "route_steps.observations.jsonl.gz",
                "routes.source.jsonl.gz",
            } or {k: v for k, v in bundle_report.items() if k != "reused"} != {
                k: v for k, v in item.items() if k != "reused"
            }:
                raise ValueError(f"Invalid intermediate bundle manifest: {bundle}")
            outputs = {bundle / name: digest for name, digest in hashes.items()}
        else:
            raise ValueError(f"Unknown intermediate category: {entry['category']}")
        for path, digest in outputs.items():
            path = path.resolve()
            if (
                path in verified or not path.is_relative_to(root)
                or not path.is_file()
                or _sha256(path) != digest
            ):
                raise ValueError(f"Invalid intermediate artifact: {path}")
            if path.name.endswith(".observations.jsonl.gz"):
                expected_inputs.add(path.resolve())
            verified[path] = digest
    if expected_inputs != {p.resolve() for p in root.rglob("*.observations.jsonl.gz")}:
        raise ValueError("Physical observation files do not match migration inventory")
    return verified


def validate_intermediate_manifest(output_root: str | Path) -> dict[str, Any]:
    """Verify all prepared output hashes against the completed migration manifest."""
    verified = verified_intermediate_artifacts(output_root)
    return {"verified_artifact_count": len(verified), "coverage_complete": True}
