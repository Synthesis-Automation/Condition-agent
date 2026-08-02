"""Deterministic one-source-file to one-intermediate-file preprocessing."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Mapping, Optional

from .models import INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
from .registry import detect_adapter, get_adapter


PREPROCESSOR_DEFINITION_VERSION = "source_preprocessor.v1.0"


@dataclass(frozen=True)
class PreprocessingProgress:
    """One progress event from file-level preprocessing."""

    phase: str
    source_file: str
    file_number: int
    file_count: int
    row_count: int
    message: str


class PreprocessingCancelled(RuntimeError):
    """Raised after an in-progress temporary artifact is safely discarded."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _output_path(source: Path, output_dir: Path) -> Path:
    return output_dir / f"{source.name}.observations.jsonl.gz"


def _cached_report(
    *,
    output_path: Path,
    source: Path,
    source_sha256: str,
    adapter_id: str,
    adapter_version: str,
) -> Optional[Dict[str, Any]]:
    """Validate and summarize an existing artifact without a sidecar log."""
    if not output_path.is_file():
        return None
    warning_counts: Counter[str] = Counter()
    kind_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    row_count = 0
    try:
        with gzip.open(output_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                record = json.loads(line)
                if not isinstance(record, Mapping):
                    return None
                provenance = record.get("source", {})
                if not isinstance(provenance, Mapping):
                    return None
                if (
                    record.get("schema_version")
                    != INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
                    or provenance.get("source_file") != source.name
                    or provenance.get("source_file_sha256") != source_sha256
                    or provenance.get("adapter_id") != adapter_id
                    or provenance.get("adapter_version") != adapter_version
                ):
                    return None
                warnings = set(record.get("warnings") or ())
                conditions = record.get("conditions") or {}
                if not isinstance(conditions, Mapping):
                    return None
                warnings.update(conditions.get("warnings") or ())
                warning_counts.update(warnings)
                kind_counts[str(record.get("observation_kind") or "")] += 1
                status_counts[str(record.get("ingestion_status") or "")] += 1
                row_count += 1
    except (EOFError, OSError, UnicodeError, json.JSONDecodeError, TypeError):
        return None
    if row_count == 0:
        return None
    return {
        "artifact_type": "source_preprocessing",
        "source_path": str(source.resolve()),
        "source_file": source.name,
        "source_sha256": source_sha256,
        "adapter_id": adapter_id,
        "adapter_version": adapter_version,
        "intermediate_schema_version": INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
        "preprocessor_definition_version": PREPROCESSOR_DEFINITION_VERSION,
        "input_row_count": row_count,
        "output_row_count": row_count,
        "observation_kind_counts": dict(sorted(kind_counts.items())),
        "ingestion_status_counts": dict(sorted(status_counts.items())),
        "warning_counts": dict(sorted(warning_counts.items())),
        "output_path": str(output_path.resolve()),
        "output_size_bytes": output_path.stat().st_size,
        "output_sha256": _sha256(output_path),
        "coverage_complete": True,
        "reused": True,
    }


def preprocess_file(
    source_path: str | Path,
    output_dir: str | Path,
    *,
    adapter_id: str = "auto",
    force: bool = False,
    progress_callback: Optional[Callable[[PreprocessingProgress], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
    file_number: int = 1,
    file_count: int = 1,
) -> Dict[str, Any]:
    """Normalize one source CSV into one deterministic compressed JSONL file."""
    source = Path(source_path)
    if not source.is_file():
        raise ValueError(f"Source file does not exist: {source}")
    if source.suffix.casefold() != ".csv":
        raise ValueError(f"Only CSV source files are supported: {source}")
    destination = Path(output_dir)
    if destination.resolve() == source.parent.resolve():
        raise ValueError("Output folder must be separate from the raw source folder")
    destination.mkdir(parents=True, exist_ok=True)
    adapter = (
        detect_adapter(source) if adapter_id == "auto" else get_adapter(adapter_id)
    )
    source_sha256 = _sha256(source)
    output_path = _output_path(source, destination)
    if not force:
        cached = _cached_report(
            output_path=output_path,
            source=source,
            source_sha256=source_sha256,
            adapter_id=adapter.adapter_id,
            adapter_version=adapter.adapter_version,
        )
        if cached is not None:
            if progress_callback is not None:
                progress_callback(
                    PreprocessingProgress(
                        "reused",
                        str(source),
                        file_number,
                        file_count,
                        int(cached["output_row_count"]),
                        f"Reused {source.name}; source and adapter are unchanged.",
                    )
                )
            return cached
    if progress_callback is not None:
        progress_callback(
            PreprocessingProgress(
                "started",
                str(source),
                file_number,
                file_count,
                0,
                f"Preprocessing {source.name} with {adapter.adapter_id}…",
            )
        )
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    warning_counts: Counter[str] = Counter()
    kind_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    row_count = 0
    try:
        with temporary.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                compresslevel=6,
                mtime=0,
            ) as compressed:
                with io.TextIOWrapper(
                    compressed, encoding="utf-8", newline="\n"
                ) as text_handle:
                    for observation in adapter.iter_observations(
                        source, source_sha256=source_sha256
                    ):
                        if cancel_check is not None and cancel_check():
                            raise PreprocessingCancelled(
                                f"Preprocessing cancelled during {source.name}"
                            )
                        payload = observation.to_dict()
                        text_handle.write(
                            json.dumps(
                                payload,
                                ensure_ascii=False,
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                        )
                        text_handle.write("\n")
                        row_count += 1
                        warning_counts.update(
                            set(observation.warnings)
                            | set(observation.conditions.warnings)
                        )
                        kind_counts[observation.observation_kind] += 1
                        status_counts[observation.ingestion_status] += 1
                        if progress_callback is not None and row_count % 1_000 == 0:
                            progress_callback(
                                PreprocessingProgress(
                                    "rows",
                                    str(source),
                                    file_number,
                                    file_count,
                                    row_count,
                                    f"{source.name}: {row_count} row(s) normalized.",
                                )
                            )
        temporary.replace(output_path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    report: Dict[str, Any] = {
        "artifact_type": "source_preprocessing",
        "source_path": str(source.resolve()),
        "source_file": source.name,
        "source_sha256": source_sha256,
        "adapter_id": adapter.adapter_id,
        "adapter_version": adapter.adapter_version,
        "intermediate_schema_version": INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
        "preprocessor_definition_version": PREPROCESSOR_DEFINITION_VERSION,
        "input_row_count": row_count,
        "output_row_count": row_count,
        "observation_kind_counts": dict(sorted(kind_counts.items())),
        "ingestion_status_counts": dict(sorted(status_counts.items())),
        "warning_counts": dict(sorted(warning_counts.items())),
        "output_path": str(output_path.resolve()),
        "output_size_bytes": output_path.stat().st_size,
        "output_sha256": _sha256(output_path),
        "coverage_complete": True,
        "reused": False,
    }
    if progress_callback is not None:
        progress_callback(
            PreprocessingProgress(
                "completed",
                str(source),
                file_number,
                file_count,
                row_count,
                f"Completed {source.name}: {row_count} observation(s).",
            )
        )
    return report


def preprocess_files(
    source_paths: Iterable[str | Path],
    output_dir: str | Path,
    *,
    adapter_id: str = "auto",
    force: bool = False,
    progress_callback: Optional[Callable[[PreprocessingProgress], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Preprocess selected files independently and return a batch summary."""
    sources = tuple(Path(path) for path in source_paths)
    if not sources:
        raise ValueError("At least one source file is required")
    duplicate_names = [
        name
        for name, count in Counter(path.name.casefold() for path in sources).items()
        if count > 1
    ]
    if duplicate_names:
        raise ValueError(
            "Selected source files have colliding names: "
            + ", ".join(sorted(duplicate_names))
        )
    reports = []
    for file_number, source in enumerate(sources, start=1):
        if cancel_check is not None and cancel_check():
            raise PreprocessingCancelled("Preprocessing cancelled between files")
        reports.append(
            preprocess_file(
                source,
                output_dir,
                adapter_id=adapter_id,
                force=force,
                progress_callback=progress_callback,
                cancel_check=cancel_check,
                file_number=file_number,
                file_count=len(sources),
            )
        )
    return {
        "artifact_type": "source_preprocessing_batch",
        "preprocessor_definition_version": PREPROCESSOR_DEFINITION_VERSION,
        "intermediate_schema_version": INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
        "output_dir": str(Path(output_dir).resolve()),
        "file_count": len(reports),
        "row_count": sum(int(report["output_row_count"]) for report in reports),
        "reused_file_count": sum(bool(report.get("reused")) for report in reports),
        "files": reports,
    }


__all__ = [
    "PREPROCESSOR_DEFINITION_VERSION",
    "PreprocessingCancelled",
    "PreprocessingProgress",
    "preprocess_file",
    "preprocess_files",
]
