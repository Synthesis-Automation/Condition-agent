"""Resumable shard compilation and deterministic full-scale library merging."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import sqlite3
import time
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass, replace
from itertools import islice
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, Literal, Sequence

from retrosynthesis_poc.chemistry import digest
from .generic_library import (
    build_generic_library,
    load_generic_library,
    save_generic_library,
    select_context_representatives,
)
from .generic_models import (
    GenericCoreTemplate,
    GenericGraphOperator,
    GenericHandleCompletionGroup,
    GenericTemplateLibrary,
)
from .retrieval_index import build_generic_retrieval_index
from .sources import iter_library_rows, source_shard_files


ProgressCallback = Callable[[Dict[str, Any]], None]


@dataclass(frozen=True)
class FullScaleBuildConfig:
    """Versioned settings that determine one reproducible shard build."""

    levels: tuple[Literal["L0", "L1", "L2"], ...] = ("L0", "L1", "L2")
    admission_mode: Literal["supported", "data_driven"] = "data_driven"
    max_precedents_per_template: int = 8
    max_rows_per_shard: int | None = None
    definition_id: str = "full_scale_operator_build.v2"

    @property
    def config_id(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return digest("FSB2", payload)


def _source_files(source: str | Path) -> tuple[Path, ...]:
    return tuple(
        sorted(
            source_shard_files(source),
            key=lambda path: (
                hashlib.sha256(path.name.encode()).hexdigest(),
                path.name,
            ),
        )
    )


def _shard_stem(source: Path) -> str:
    readable = "".join(
        character if character.isalnum() else "_"
        for character in source.name[:48]
    ).strip("_")
    identity = hashlib.sha256(str(source.resolve()).encode()).hexdigest()[:16]
    return f"{readable}.{identity}"


def _source_identity(source: Path) -> Dict[str, Any]:
    stat = source.stat()
    return {
        "source_path": str(source.resolve()),
        "source_size": stat.st_size,
        "source_mtime_ns": stat.st_mtime_ns,
    }


def _iter_shard_rows(
    source: Path,
    max_rows: int | None,
) -> Iterator[Dict[str, Any]]:
    rows: Iterable[Dict[str, Any]] = iter_library_rows(source)
    if max_rows is not None:
        rows = islice(rows, max_rows)
    for ordinal, row in enumerate(rows, start=1):
        yield {
            **row,
            "_build_source_shard": source.name,
            "_build_source_row_number": ordinal,
        }


def _write_json(path: Path, value: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _manifest_current(
    manifest_path: Path,
    library_path: Path,
    ledger_path: Path,
    expected: Dict[str, Any],
) -> bool:
    if not manifest_path.is_file() or not library_path.is_file() or not ledger_path.is_file():
        return False
    try:
        value = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return all(value.get(key) == item for key, item in expected.items())


def compile_operator_shard(
    source: str | Path,
    output_directory: str | Path,
    *,
    config: FullScaleBuildConfig | None = None,
    force: bool = False,
) -> Dict[str, Any]:
    """Compile one shard to resumable library, ledger, and manifest artifacts."""

    settings = config or FullScaleBuildConfig()
    source_path = Path(source)
    output = Path(output_directory)
    shard_directory = output / "shards"
    shard_directory.mkdir(parents=True, exist_ok=True)
    stem = _shard_stem(source_path)
    library_path = shard_directory / f"{stem}.library.json.gz"
    ledger_path = shard_directory / f"{stem}.admissions.jsonl.gz"
    manifest_path = shard_directory / f"{stem}.manifest.json"
    expected = {
        **_source_identity(source_path),
        "config_id": settings.config_id,
    }
    if not force and _manifest_current(
        manifest_path,
        library_path,
        ledger_path,
        expected,
    ):
        return {
            **json.loads(manifest_path.read_text(encoding="utf-8")),
            "_reused": True,
        }

    temporary_ledger = ledger_path.with_suffix(ledger_path.suffix + ".tmp")
    admission_counts: Counter[str] = Counter()
    with temporary_ledger.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8") as text_handle:

                def record_admission(record: Dict[str, Any]) -> None:
                    admission_counts[
                        str(record.get("reason") or record.get("status") or "unknown")
                    ] += 1
                    text_handle.write(
                        json.dumps(
                            record,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        + "\n"
                    )

                library = build_generic_library(
                    _iter_shard_rows(
                        source_path,
                        settings.max_rows_per_shard,
                    ),
                    levels=settings.levels,
                    admission_mode=settings.admission_mode,
                    max_precedents_per_template=(
                        settings.max_precedents_per_template
                    ),
                    admission_callback=record_admission,
                )

    temporary_library = library_path.with_suffix(library_path.suffix + ".tmp")
    save_generic_library(library, temporary_library)
    temporary_ledger.replace(ledger_path)
    temporary_library.replace(library_path)
    manifest = {
        **expected,
        "definition_id": settings.definition_id,
        "source_rows": library.source_row_count,
        "accepted_observations": library.accepted_observation_count,
        "template_count": len(library.templates),
        "operator_count": len(library.operators),
        "completion_group_count": len(library.completion_groups),
        "admission_counts": dict(sorted(admission_counts.items())),
        "library_path": str(library_path.resolve()),
        "ledger_path": str(ledger_path.resolve()),
    }
    _write_json(manifest_path, manifest)
    return {**manifest, "_reused": False}


def _ledger_rows(path: Path) -> Iterator[Dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            value = json.loads(line)
            if isinstance(value, dict):
                yield value


def _initialize_support_database(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        DROP TABLE IF EXISTS accepted;
        DROP TABLE IF EXISTS template_support;
        DROP TABLE IF EXISTS operator_support;
        DROP TABLE IF EXISTS completion_support;
        CREATE TABLE accepted (
            observation_key TEXT PRIMARY KEY,
            reference_key TEXT NOT NULL,
            annotated INTEGER NOT NULL
        );
        CREATE TABLE template_support (
            template_id TEXT NOT NULL,
            observation_key TEXT NOT NULL,
            reference_key TEXT NOT NULL,
            PRIMARY KEY (template_id, observation_key)
        );
        CREATE TABLE operator_support (
            operator_id TEXT NOT NULL,
            observation_key TEXT NOT NULL,
            reference_key TEXT NOT NULL,
            PRIMARY KEY (operator_id, observation_key)
        );
        CREATE TABLE completion_support (
            completion_group_id TEXT NOT NULL,
            observation_key TEXT NOT NULL,
            reference_key TEXT NOT NULL,
            PRIMARY KEY (completion_group_id, observation_key)
        );
        """
    )


def _insert_admission_support(
    connection: sqlite3.Connection,
    record: Dict[str, Any],
) -> bool:
    observation_key = str(
        record.get("observation_key") or record.get("reaction_id") or ""
    )
    if not observation_key:
        return False
    reference_key = str(record.get("support_key") or observation_key)
    cursor = connection.execute(
        "INSERT OR IGNORE INTO accepted VALUES (?, ?, ?)",
        (
            observation_key,
            reference_key,
            int(bool(record.get("named_annotations"))),
        ),
    )
    if cursor.rowcount == 0:
        return False
    connection.executemany(
        "INSERT OR IGNORE INTO template_support VALUES (?, ?, ?)",
        (
            (str(template_id), observation_key, reference_key)
            for template_id in record.get("template_ids") or ()
        ),
    )
    connection.executemany(
        "INSERT OR IGNORE INTO operator_support VALUES (?, ?, ?)",
        (
            (str(operator_id), observation_key, reference_key)
            for operator_id in record.get("operator_ids") or ()
        ),
    )
    connection.executemany(
        "INSERT OR IGNORE INTO completion_support VALUES (?, ?, ?)",
        (
            (str(group_id), observation_key, reference_key)
            for group_id in record.get("completion_group_ids") or ()
        ),
    )
    return True


def _support_counts(
    connection: sqlite3.Connection,
    table: str,
    identity_column: str,
) -> dict[str, tuple[int, int]]:
    query = (
        f"SELECT {identity_column}, COUNT(*), COUNT(DISTINCT reference_key) "
        f"FROM {table} GROUP BY {identity_column}"
    )
    return {
        str(identity): (int(observations), int(references))
        for identity, observations, references in connection.execute(query)
    }


def _merge_template(
    current: GenericCoreTemplate,
    incoming: GenericCoreTemplate,
    max_precedents: int,
) -> GenericCoreTemplate:
    return replace(
        current,
        precedents=select_context_representatives(
            (*current.precedents, *incoming.precedents),
            max_precedents,
        ),
    )


def merge_operator_shards(
    manifests: Sequence[Dict[str, Any]],
    output_directory: str | Path,
    *,
    config: FullScaleBuildConfig | None = None,
    progress_callback: ProgressCallback | None = None,
    progress_started_at: float | None = None,
) -> tuple[GenericTemplateLibrary, Dict[str, Any]]:
    """Merge shard artifacts with exact observation/reference deduplication."""

    settings = config or FullScaleBuildConfig()
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    database_path = output / "support.sqlite3"
    connection = sqlite3.connect(database_path)
    _initialize_support_database(connection)
    templates: dict[str, GenericCoreTemplate] = {}
    rejection_counts: Counter[str] = Counter()
    source_rows = 0
    duplicate_accepted = 0
    merge_started = progress_started_at or time.monotonic()
    ordered_manifests = sorted(
        manifests,
        key=lambda value: value["source_path"],
    )
    merge_progress_step = max(1, len(ordered_manifests) // 20)
    if progress_callback is not None:
        progress_callback(
            {
                "phase": "merge",
                "completed_shards": 0,
                "total_shards": len(ordered_manifests),
                "source_rows": 0,
                "elapsed_seconds": time.monotonic() - merge_started,
            }
        )
    try:
        for ordinal, manifest in enumerate(ordered_manifests, start=1):
            source_rows += int(manifest.get("source_rows") or 0)
            partial = load_generic_library(str(manifest["library_path"]))
            for template in partial.templates:
                current = templates.get(template.template_id)
                templates[template.template_id] = (
                    template
                    if current is None
                    else _merge_template(
                        current,
                        template,
                        settings.max_precedents_per_template,
                    )
                )
            for record in _ledger_rows(Path(str(manifest["ledger_path"]))):
                if record.get("status") == "accepted":
                    if not _insert_admission_support(connection, record):
                        duplicate_accepted += 1
                else:
                    rejection_counts[
                        str(record.get("reason") or "unknown_rejection")
                    ] += 1
            if progress_callback is not None and (
                ordinal % merge_progress_step == 0
                or ordinal == len(ordered_manifests)
            ):
                progress_callback(
                    {
                        "phase": "merge",
                        "completed_shards": ordinal,
                        "total_shards": len(ordered_manifests),
                        "source_rows": source_rows,
                        "template_count": len(templates),
                        "elapsed_seconds": time.monotonic() - merge_started,
                    }
                )
        connection.commit()
        if progress_callback is not None:
            progress_callback(
                {
                    "phase": "finalize",
                    "completed_shards": len(ordered_manifests),
                    "total_shards": len(ordered_manifests),
                    "source_rows": source_rows,
                    "template_count": len(templates),
                    "elapsed_seconds": time.monotonic() - merge_started,
                }
            )
        template_counts = _support_counts(
            connection,
            "template_support",
            "template_id",
        )
        operator_counts = _support_counts(
            connection,
            "operator_support",
            "operator_id",
        )
        completion_counts = _support_counts(
            connection,
            "completion_support",
            "completion_group_id",
        )
        accepted_count = int(
            connection.execute("SELECT COUNT(*) FROM accepted").fetchone()[0]
        )
        annotated_count = int(
            connection.execute(
                "SELECT COUNT(*) FROM accepted WHERE annotated = 1"
            ).fetchone()[0]
        )
    finally:
        connection.close()

    finalized = tuple(
        replace(
            template,
            observation_support=template_counts.get(template.template_id, (0, 0))[0],
            independent_reference_support=template_counts.get(
                template.template_id,
                (0, 0),
            )[1],
        )
        for template in sorted(templates.values(), key=lambda item: item.template_id)
        if template.template_id in template_counts
    )
    by_operator: dict[str, list[GenericCoreTemplate]] = {}
    by_completion: dict[str, list[GenericCoreTemplate]] = {}
    for template in finalized:
        by_operator.setdefault(template.operator_id, []).append(template)
        completion_id = digest(
            "COMP2",
            template.operator_id,
            template.handle_signature,
        )
        by_completion.setdefault(completion_id, []).append(template)
    operators = tuple(
        GenericGraphOperator(
            operator_id=operator_id,
            operator_signature=members[0].operator_signature,
            edit_tokens=tuple(
                sorted({token for member in members for token in member.edit_tokens})
            ),
            realization_ids=tuple(
                sorted({member.realization_id for member in members})
            ),
            abstraction_levels=tuple(
                sorted({member.abstraction_level for member in members})
            ),
            observation_support=operator_counts.get(operator_id, (0, 0))[0],
            independent_reference_support=operator_counts.get(operator_id, (0, 0))[1],
            named_annotations=tuple(
                sorted(
                    {
                        annotation
                        for member in members
                        for annotation in member.named_annotations
                    }
                )
            ),
        )
        for operator_id, members in sorted(by_operator.items())
    )
    completion_groups = tuple(
        GenericHandleCompletionGroup(
            completion_group_id=completion_id,
            operator_id=members[0].operator_id,
            completion_signature=members[0].handle_signature,
            synthon_signatures=tuple(
                sorted({member.synthon_signature for member in members})
            ),
            realization_ids=tuple(
                sorted({member.realization_id for member in members})
            ),
            template_ids=tuple(sorted(member.template_id for member in members)),
            handle_signatures=tuple(
                sorted({member.handle_signature for member in members})
            ),
            observation_support=completion_counts.get(completion_id, (0, 0))[0],
            independent_reference_support=completion_counts.get(
                completion_id,
                (0, 0),
            )[1],
        )
        for completion_id, members in sorted(by_completion.items())
    )
    library = GenericTemplateLibrary(
        templates=finalized,
        source_row_count=source_rows,
        accepted_observation_count=accepted_count,
        rejection_counts=dict(sorted(rejection_counts.items())),
        definition={
            "definition_id": "generic_core_retrosynthesis_full_scale.v3",
            "build_config_id": settings.config_id,
            "routing_source": "normalized_graph_edits",
            "admission_mode": settings.admission_mode,
            "levels": list(settings.levels),
            "source_shard_count": len(manifests),
            "annotated_accepted_observation_count": annotated_count,
            "unannotated_accepted_observation_count": (
                accepted_count - annotated_count
            ),
            "deduplicated_accepted_observation_count": accepted_count,
            "duplicate_accepted_observation_count": duplicate_accepted,
            "context_representative_policy": "deterministic_distinct_context_bins",
        },
        operators=operators,
        completion_groups=completion_groups,
    )
    library = replace(
        library,
        retrieval_index=build_generic_retrieval_index(library),
    )
    library_path = output / "operator_library_v3.json.gz"
    save_generic_library(library, library_path)
    report = {
        "definition_id": "full_scale_operator_build_report.v2",
        "build_config_id": settings.config_id,
        "source_shards": len(manifests),
        "source_rows": source_rows,
        "accepted_observations": accepted_count,
        "acceptance_fraction": accepted_count / source_rows if source_rows else 0.0,
        "duplicate_accepted_observations": duplicate_accepted,
        "template_count": len(finalized),
        "operator_count": len(operators),
        "completion_group_count": len(completion_groups),
        "realization_count": len(
            {template.realization_id for template in finalized}
        ),
        "rejection_counts": dict(sorted(rejection_counts.items())),
        "library_path": str(library_path.resolve()),
        "support_database_path": str(database_path.resolve()),
    }
    _write_json(output / "build_report.json", report)
    if progress_callback is not None:
        progress_callback(
            {
                "phase": "complete",
                **report,
                "completed_shards": len(ordered_manifests),
                "total_shards": len(ordered_manifests),
                "elapsed_seconds": time.monotonic() - merge_started,
            }
        )
    return library, report


def build_full_scale_operator_library(
    source: str | Path,
    output_directory: str | Path,
    *,
    config: FullScaleBuildConfig | None = None,
    max_shards: int | None = None,
    workers: int = 1,
    force: bool = False,
    progress_callback: ProgressCallback | None = None,
    progress_interval_seconds: float = 30.0,
) -> tuple[GenericTemplateLibrary, Dict[str, Any]]:
    """Compile selected shards resumably and merge the resulting v3 library."""

    settings = config or FullScaleBuildConfig()
    files = _source_files(source)
    if max_shards is not None:
        if max_shards < 1:
            raise ValueError("max shards must be positive")
        files = files[:max_shards]
    if not files:
        raise ValueError("source contains no JSONL gzip shards")
    if workers < 1:
        raise ValueError("workers must be positive")
    if progress_interval_seconds <= 0:
        raise ValueError("progress interval must be positive")
    started = time.monotonic()
    completed_manifests: list[Dict[str, Any]] = []

    def emit_compile_progress(
        *,
        newly_completed_shards: int = 0,
        newly_reused_shards: int = 0,
    ) -> None:
        if progress_callback is None:
            return
        source_rows = sum(
            int(manifest.get("source_rows") or 0)
            for manifest in completed_manifests
        )
        accepted = sum(
            int(manifest.get("accepted_observations") or 0)
            for manifest in completed_manifests
        )
        progress_callback(
            {
                "phase": "compile",
                "completed_shards": len(completed_manifests),
                "total_shards": len(files),
                "source_rows": source_rows,
                "accepted_observations": accepted,
                "reused_shards": sum(
                    int(bool(manifest.get("_reused")))
                    for manifest in completed_manifests
                ),
                "newly_completed_shards": newly_completed_shards,
                "newly_reused_shards": newly_reused_shards,
                "active_shards": min(
                    workers,
                    max(0, len(files) - len(completed_manifests)),
                ),
                "queued_shards": max(
                    0,
                    len(files) - len(completed_manifests) - workers,
                ),
                "workers": workers,
                "elapsed_seconds": time.monotonic() - started,
            }
        )

    emit_compile_progress()
    if workers == 1:
        last_progress = started
        newly_completed = 0
        newly_reused = 0
        for ordinal, path in enumerate(files, start=1):
            manifest = compile_operator_shard(
                path,
                output_directory,
                config=settings,
                force=force,
            )
            completed_manifests.append(manifest)
            newly_completed += 1
            newly_reused += int(bool(manifest.get("_reused")))
            now = time.monotonic()
            if (
                now - last_progress >= progress_interval_seconds
                or ordinal == len(files)
            ):
                emit_compile_progress(
                    newly_completed_shards=newly_completed,
                    newly_reused_shards=newly_reused,
                )
                last_progress = now
                newly_completed = 0
                newly_reused = 0
    else:
        jobs = {
            path: (path, Path(output_directory), settings, force)
            for path in files
        }
        with ProcessPoolExecutor(max_workers=workers) as executor:
            pending = {
                executor.submit(_compile_operator_shard_job, job): path
                for path, job in jobs.items()
            }
            last_progress = time.monotonic()
            newly_completed = 0
            newly_reused = 0
            while pending:
                done, not_done = wait(
                    pending,
                    timeout=progress_interval_seconds,
                    return_when=FIRST_COMPLETED,
                )
                pending = set(not_done)
                completed_batch = []
                for future in done:
                    completed_batch.append(future.result())
                completed_manifests.extend(completed_batch)
                newly_completed += len(completed_batch)
                newly_reused += sum(
                    int(bool(manifest.get("_reused")))
                    for manifest in completed_batch
                )
                now = time.monotonic()
                if (
                    progress_callback is not None
                    and (
                        now - last_progress >= progress_interval_seconds
                        or not pending
                    )
                ):
                    emit_compile_progress(
                        newly_completed_shards=newly_completed,
                        newly_reused_shards=newly_reused,
                    )
                    last_progress = now
                    newly_completed = 0
                    newly_reused = 0
    manifests = tuple(completed_manifests)
    return merge_operator_shards(
        manifests,
        output_directory,
        config=settings,
        progress_callback=progress_callback,
        progress_started_at=started,
    )


def _compile_operator_shard_job(
    job: tuple[Path, Path, FullScaleBuildConfig, bool],
) -> Dict[str, Any]:
    source, output, config, force = job
    return compile_operator_shard(
        source,
        output,
        config=config,
        force=force,
    )


__all__ = [
    "FullScaleBuildConfig",
    "build_full_scale_operator_library",
    "compile_operator_shard",
    "merge_operator_shards",
    "ProgressCallback",
]
