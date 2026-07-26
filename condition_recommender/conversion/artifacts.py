"""Build scalable recommendation artifacts and their human review view."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, Mapping, Optional

from ..generic_indexing import build_generic_index, save_generic_index
from .concise_review import (
    ConciseReviewConversionCancelled,
    ConciseReviewProgress,
    export_concise_reaction_review_csv,
    iter_canonical_records,
)
from .sharded import (
    ShardedConversionCancelled,
    ShardedConversionProgress,
    convert_datasets_sharded,
)

RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION = "1.0"


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


def _validate_paths(dataset_path: Path, output_dir: Path) -> None:
    source = dataset_path.resolve()
    destination = output_dir.resolve()
    if source.is_dir() and (destination == source or source in destination.parents):
        raise ValueError(
            "Output folder must be outside the source dataset folder so "
            "generated CSV files are not discovered as source data."
        )


def build_recommendation_artifacts(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    shard_size: int = 1_000,
    workers: int = 1,
    build_fast_index: bool = True,
    checkpoint_interval: int = 1,
    progress_callback: Optional[
        Callable[[RecommendationArtifactProgress], None]
    ] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Build restartable canonical data, a review CSV, and a fast index.

    Canonical records remain the source of truth. The persisted generic index
    contains only retrieval-eligible fields and lookup maps, while the CSV is a
    concise review projection and is never used as recommendation input.
    """
    source = Path(dataset_path)
    destination = Path(output_dir)
    _validate_paths(source, destination)
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

    try:
        conversion_report = convert_datasets_sharded(
            source,
            destination,
            shard_size=shard_size,
            mode="full",
            workers=workers,
            checkpoint_interval=checkpoint_interval,
            merge_records=False,
            progress_callback=on_conversion_progress,
            cancel_check=cancel_check,
        )
    except ShardedConversionCancelled as exc:
        raise RecommendationArtifactBuildCancelled(str(exc)) from exc

    failed_shards = int(conversion_report.get("failed_shard_count") or 0)
    if failed_shards:
        raise RuntimeError(
            f"Canonical conversion has {failed_shards} failed shard(s); "
            "review and index artifacts were not built."
        )
    records_path = destination / "shard_manifest.json"
    total_rows = int(conversion_report.get("output_row_count") or 0)
    if cancel_check is not None and cancel_check():
        raise RecommendationArtifactBuildCancelled(
            "Build cancelled after canonical data completed; run again with "
            "the same output folder to reuse it."
        )

    review_path = destination / "reaction_review.csv"
    notify(
        "review_started",
        "Writing the concise chemistry review CSV…",
        row_count=0,
    )

    def on_review_progress(progress: ConciseReviewProgress) -> None:
        notify(
            "review_completed"
            if progress.phase == "completed"
            else "review_rows",
            progress.message,
            row_count=progress.row_count,
    )

    try:
        export_concise_reaction_review_csv(
            records_path,
            review_path,
            progress_callback=on_review_progress,
            cancel_check=cancel_check,
        )
    except ConciseReviewConversionCancelled as exc:
        raise RecommendationArtifactBuildCancelled(
            "Build cancelled during review export; canonical shards remain "
            "available for resume."
        ) from exc

    index_report: Optional[Dict[str, Any]] = None
    index_path = destination / "generic_index.json.gz"
    if build_fast_index:
        if cancel_check is not None and cancel_check():
            raise RecommendationArtifactBuildCancelled(
                "Build cancelled before index creation; canonical data and "
                "the review CSV are complete."
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
                        "and the review CSV are complete."
                    )
                scanned_rows += 1
                if scanned_rows % 1_000 == 0:
                    notify(
                        "index_rows",
                        f"Indexed input scan: {scanned_rows} record(s).",
                        row_count=scanned_rows,
                    )
                yield record

        index = build_generic_index(records())
        if cancel_check is not None and cancel_check():
            raise RecommendationArtifactBuildCancelled(
                "Build cancelled before saving the fast-load index."
            )
        index_report = save_generic_index(index, index_path)
        notify(
            "index_completed",
            (
                "Fast-load index complete: "
                f"{index_report['row_count']} eligible precedent(s)."
            ),
            row_count=int(index_report["row_count"]),
        )

    artifacts: Dict[str, Any] = {
        "canonical_manifest": _artifact_entry(records_path, destination),
        "review_csv": _artifact_entry(review_path, destination),
    }
    if index_report is not None:
        artifacts["fast_index"] = _artifact_entry(index_path, destination)
    shard_paths = tuple((destination / "shards").glob("*.jsonl.gz"))
    shard_size_bytes = sum(path.stat().st_size for path in shard_paths)
    all_files = tuple(path for path in destination.rglob("*") if path.is_file())
    warnings = []
    if not build_fast_index and index_path.is_file():
        warnings.append(
            "An older generic_index.json.gz exists but was not rebuilt; do "
            "not use it unless it matches the current canonical records."
        )
    report: Dict[str, Any] = {
        "schema_version": RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION,
        "artifact_type": "recommendation_artifact_build",
        "source_path": str(source.resolve()),
        "output_dir": str(destination.resolve()),
        "settings": {
            "shard_size": shard_size,
            "workers": workers,
            "build_fast_index": build_fast_index,
            "compression": "gzip",
        },
        "source_file_count": int(conversion_report["source_file_count"]),
        "record_count": total_rows,
        "eligible_index_record_count": (
            int(index_report["row_count"]) if index_report is not None else None
        ),
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
    "RECOMMENDATION_ARTIFACT_WORKFLOW_SCHEMA_VERSION",
    "RecommendationArtifactBuildCancelled",
    "RecommendationArtifactProgress",
    "build_recommendation_artifacts",
]
