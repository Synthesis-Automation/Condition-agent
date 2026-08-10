"""Resolve canonical Full and Compact recommendation-library inputs."""

from __future__ import annotations

import fnmatch
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Literal

from retrosynthesis_poc.library import iter_rows


LibraryMode = Literal["full", "compact"]
LIBRARY_MODES: tuple[LibraryMode, ...] = ("full", "compact")
COMBINED_RECORDS_FILENAME = "combined_records.jsonl.gz"
COMBINED_BATCH_MANIFEST_FILENAME = "combined_batch_manifest.json"


def resolve_library_mode(
    source: str | Path,
    mode: str = "full",
) -> Path:
    """Resolve a recommendation-library root to one mode-specific directory.

    Explicit files, legacy shard directories, and already selected ``full`` or
    ``compact`` directories remain valid inputs.
    """

    normalized = mode.strip().casefold()
    if normalized not in LIBRARY_MODES:
        raise ValueError(f"unsupported library mode: {mode}")
    root = Path(source)
    if root.is_file() or root.name.casefold() in LIBRARY_MODES:
        return root
    candidate = root / normalized
    if candidate.is_dir():
        return candidate
    return root


def combined_records_source(source: str | Path) -> Path:
    """Return the deduplicated corpus for a mode directory when available."""

    root = Path(source)
    candidate = root / COMBINED_RECORDS_FILENAME
    return candidate if root.is_dir() and candidate.is_file() else root


def _matches_include(row: Dict[str, Any], patterns: tuple[str, ...]) -> bool:
    if not patterns:
        return True
    source_name = Path(str(row.get("source_path") or "")).name.casefold()
    return any(
        fnmatch.fnmatchcase(source_name, pattern.casefold())
        for pattern in patterns
    )


def iter_library_rows(
    source: str | Path,
    *,
    include: Iterable[str] = (),
) -> Iterator[Dict[str, Any]]:
    """Yield canonical rows without duplicating batch and combined artifacts."""

    selected = combined_records_source(source)
    patterns = tuple(include)
    if selected.is_file() and selected.name == COMBINED_RECORDS_FILENAME:
        for row in iter_rows(selected):
            if _matches_include(row, patterns):
                yield row
        return
    yield from iter_rows(selected, include=patterns)


def _manifest_path(mode_root: Path, configured_path: object) -> Path:
    candidate = Path(str(configured_path or ""))
    if candidate.is_file():
        return candidate
    if candidate.parent.name:
        relocated = (
            mode_root
            / "batches"
            / candidate.parent.name
            / "shard_manifest.json"
        )
        if relocated.is_file():
            return relocated
    raise FileNotFoundError(f"saved batch manifest is unavailable: {candidate}")


def _combined_batch_manifests(mode_root: Path) -> tuple[Path, ...]:
    combined_path = mode_root / COMBINED_BATCH_MANIFEST_FILENAME
    if not combined_path.is_file():
        return ()
    value = json.loads(combined_path.read_text(encoding="utf-8"))
    return tuple(
        _manifest_path(mode_root, entry.get("path"))
        for entry in value.get("batch_manifests") or ()
        if isinstance(entry, dict)
    )


def source_shard_files(source: str | Path) -> tuple[Path, ...]:
    """Return canonical physical shards represented by the selected library.

    A mode-specific combined manifest is authoritative, preventing incomplete
    or newly saved-but-not-combined batches from leaking into a build.
    """

    root = Path(source)
    if root.is_file():
        return (root,)
    manifests = _combined_batch_manifests(root)
    if manifests:
        files: list[Path] = []
        seen: set[Path] = set()
        for manifest_path in manifests:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            source_files = tuple(manifest.get("source_files") or ())
            if not source_files or any(
                not item.get("coverage_complete") for item in source_files
            ):
                raise ValueError(f"saved batch is incomplete: {manifest_path}")
            for entry in manifest.get("shards") or ():
                if entry.get("status") != "complete":
                    raise ValueError(
                        f"saved batch contains an incomplete shard: {manifest_path}"
                    )
                shard = manifest_path.parent / str(entry.get("output_path") or "")
                if not shard.is_file():
                    raise FileNotFoundError(shard)
                resolved = shard.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    files.append(shard)
        return tuple(files)
    combined = root / COMBINED_RECORDS_FILENAME
    if combined.is_file():
        return (combined,)
    return tuple(sorted(root.glob("*.jsonl.gz"), key=lambda path: path.name))


__all__ = [
    "COMBINED_BATCH_MANIFEST_FILENAME",
    "COMBINED_RECORDS_FILENAME",
    "LIBRARY_MODES",
    "LibraryMode",
    "combined_records_source",
    "iter_library_rows",
    "resolve_library_mode",
    "source_shard_files",
]
