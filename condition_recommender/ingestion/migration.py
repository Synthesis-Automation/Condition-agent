"""Rebuild a consolidated intermediate corpus from the local raw source tree."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from .artifacts import PreprocessingProgress, _sha256, preprocess_file
from .registry import detect_adapter
from .route_release import prepare_route_release

MIGRATION_DEFINITION_VERSION = "intermediate_corpus_migration.v1"


def regenerate_intermediate_datasets(
    raw_root: str | Path,
    output_root: str | Path,
    *,
    force: bool = False,
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
    routes = sorted(source.rglob("*.jsonl.gz"))
    # Preflight everything before writing a large corpus.
    adapters = {path: detect_adapter(path) for path in csvs}
    for path in routes:
        if path.name != "uspto.higher-level.routes.jsonl.gz":
            raise ValueError(f"Unsupported route source: {path}")
        if not (path.parent.parent / "reactions" / "uspto_original.csv").is_file():
            raise ValueError(f"Missing original reaction export for {path}")
    superseded = {
        path.parent.parent / "reactions" / "uspto_original.csv" for path in routes
    }
    entries = []
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
        report = preprocess_file(
            path,
            output,
            force=force,
            progress_callback=progress_callback,
        )
        entries.append(
            {
                "source_relative_path": path.relative_to(source).as_posix(),
                "category": category,
                "report": report,
            }
        )
    route_reports = []
    for path in routes:
        original = path.parent.parent / "reactions" / "uspto_original.csv"
        output = (
            destination / "routes" / path.parent.parent.relative_to(source / "routes")
        )
        report = prepare_route_release(
            path,
            original,
            output,
            force=force,
            progress_callback=progress_callback,
        )
        route_reports.append(report)
        entries.append(
            {
                "source_relative_path": path.relative_to(source).as_posix(),
                "category": "routes_and_single_step",
                "report": report,
            }
        )
        entries.append(
            {
                "source_relative_path": original.relative_to(source).as_posix(),
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
        "source_dataset_count": len(entries),
        "sources": entries,
        "single_step_observation_count": sum(
            entry["report"]["output_row_count"]
            for entry in entries
            if entry["category"] == "single_step"
        )
        + sum(r["counts"]["output_step_observations"] for r in route_reports),
        "abstraction_count": sum(
            entry["report"]["output_row_count"]
            for entry in entries
            if entry["category"] == "abstractions"
        ),
        "route_count": sum(r["counts"]["route_count"] for r in route_reports),
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
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(target)
    return report


def validate_intermediate_manifest(output_root: str | Path) -> dict[str, Any]:
    """Verify all prepared output hashes against the completed migration manifest."""
    root = Path(output_root).resolve()
    report = json.loads((root / "migration.manifest.json").read_text(encoding="utf-8"))
    verified = 0
    for entry in report["sources"]:
        item = entry.get("report")
        if item is None:
            continue
        outputs = (
            {Path(item["output_path"]): item["output_sha256"]}
            if entry["category"] in {"single_step", "abstractions"}
            else {
                root
                / "routes"
                / Path(entry["source_relative_path"]).parent.parent.relative_to(
                    "routes"
                )
                / name: digest
                for name, digest in item["output_sha256"].items()
            }
        )
        for path, digest in outputs.items():
            if not path.is_relative_to(root) or _sha256(path) != digest:
                raise ValueError(f"Invalid intermediate artifact: {path}")
            verified += 1
    if not report["coverage_complete"]:
        raise ValueError("Migration coverage is incomplete")
    return {"verified_artifact_count": verified, "coverage_complete": True}
