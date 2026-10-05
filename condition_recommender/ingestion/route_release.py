"""Lossless preparation of the released higher-level route source.

This module observes source fields only. It neither reconstructs route order nor
validates chemistry. Abstractions remain weak source annotations. Original step
IDs are deduplicated with a disk-backed ledger while all route occurrences remain
in a separate archive. The route JSONL retains agents omitted by the CSV export.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import sqlite3
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Iterator, TextIO

from .adapters.base import observation_id, supplied_mapping_status
from .artifacts import PreprocessingCancelled, PreprocessingProgress, _sha256
from .models import (
    CanonicalSourceObservation,
    ConditionComponentClaim,
    ConditionInput,
    ConditionStageInput,
    ReactionEvidenceInput,
    SourceIdentifier,
    SourceProvenance,
)

ROUTE_RELEASE_DEFINITION_VERSION = "higher_level_route_source.v1"


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _text_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


@contextmanager
def _compressed_writer(path: Path) -> Iterator[TextIO]:
    """Write reproducible gzip bytes through an atomic temporary file."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="\n") as stream:
                    yield stream
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _reaction_fields(reaction: str) -> tuple[str, str]:
    parts = reaction.split(">")
    if len(parts) != 3 or not parts[0] or not parts[2]:
        return reaction, ""
    return f"{parts[0]}>>{parts[2]}", parts[1]


def _step_observation(
    *,
    source: Path,
    source_hash: str,
    source_id: str,
    reaction: str,
    line: int,
    route_id: str,
    subtree_id: str,
    entry: dict[str, Any],
    occurrences: int,
    conflicting: bool,
) -> CanonicalSourceObservation:
    """Separate supplied agent structures without inventing their chemical roles."""
    transformation, middle = _reaction_fields(reaction)
    components = tuple(
        ConditionComponentClaim(
            component_key=f"agent_{index}",
            source_slot="reaction_middle_field",
            source_role_hint=None,
            identifiers=(SourceIdentifier("smiles", value, "reaction_middle_field"),),
            provenance={"source_declared_role": "unresolved", "source_position": index},
        )
        for index, value in enumerate(filter(None, middle.split(".")), start=1)
    )
    warnings = []
    parts = reaction.split(">")
    if len(parts) != 3 or not parts[0] or not parts[2]:
        warnings.append("MISSING_OR_INCOMPLETE_REACTION_SMILES")
    if conflicting:
        warnings.append("CONFLICTING_REACTION_STRINGS_FOR_SOURCE_ID")
    patent = route_id.rsplit("_", 1)[0]
    record_id = source_id
    if conflicting:
        record_id += f":variant-{_text_hash(reaction)}"
    return CanonicalSourceObservation(
        observation_id=observation_id(
            adapter_id=ROUTE_RELEASE_DEFINITION_VERSION,
            source_sha256=source_hash,
            row_number=line,
            record_id=record_id,
        ),
        observation_kind="structure_backed",
        source=SourceProvenance(
            corpus_id="higher_level_original_route_steps",
            release_id=source.name,
            adapter_id=ROUTE_RELEASE_DEFINITION_VERSION,
            adapter_version="1.0",
            source_file=source.name,
            source_file_sha256=source_hash,
            source_row_number=line,
            source_record_id=record_id,
            source_groups={
                "patent_id": patent,
                "first_route_id": route_id,
                "first_subtree_id": subtree_id,
                "original_reaction_id": source_id,
                "source_occurrence_count": str(occurrences),
                "route_membership_artifact": "routes.source.jsonl.gz",
            },
            reference=patent,
        ),
        reaction=ReactionEvidenceInput(
            evidence_kind="source_structure",
            reaction_smiles=transformation or None,
            supplied_mapping_status=supplied_mapping_status(transformation),
            structure_available=bool(transformation),
            source_labels={
                "source_reaction_smiles": reaction,
                "route_order": "unresolved_source_array_order_is_not_chronology",
            },
        ),
        conditions=ConditionInput(
            components=components,
            stages=(
                ConditionStageInput(
                    stage_index=1,
                    component_keys=tuple(x.component_key for x in components),
                    provenance={"source": "reaction_middle_field"},
                ),
            )
            if components
            else (),
            warnings=("REACTION_MIDDLE_FIELD_ROLES_UNRESOLVED",)
            if middle
            else ("CONDITION_DETAILS_UNAVAILABLE",),
        ),
        ingestion_status="review" if warnings else "accepted",
        warnings=tuple(warnings),
        raw_fields=dict(entry),
    )


def prepare_route_release(
    route_source: str | Path,
    original_csv: str | Path,
    output_dir: str | Path,
    *,
    progress_callback: Callable[[PreprocessingProgress], None] | None = None,
    force: bool = False,
    cancel_check: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Extract original steps and preserve every raw route; cross-check CSV coverage.

    Repeated IDs with conflicting strings retain every distinct variant for review.
    Identical strings from different source IDs remain separate observations.
    Missing CSV coverage or inconsistent CSV structures fail the migration gate.
    """
    source, csv_source, destination = (
        Path(route_source),
        Path(original_csv),
        Path(output_dir),
    )
    destination.mkdir(parents=True, exist_ok=True)
    source_hash, csv_hash = _sha256(source), _sha256(csv_source)
    manifest = destination / "route_release.manifest.json"
    if manifest.exists() and not force:
        previous = json.loads(manifest.read_text(encoding="utf-8"))
        if (
            previous.get("definition_version") == ROUTE_RELEASE_DEFINITION_VERSION
            and previous.get("source_sha256") == source_hash
            and previous.get("original_csv_sha256") == csv_hash
            and all(
                (destination / name).is_file() and _sha256(destination / name) == digest
                for name, digest in previous.get("output_sha256", {}).items()
            )
            and previous.get("coverage_complete") is True
        ):
            return {**previous, "reused": True}

    def emit(phase: str, count: int, message: str) -> None:
        if progress_callback:
            progress_callback(
                PreprocessingProgress(phase, str(source), 1, 1, count, message)
            )

    def check_cancelled() -> None:
        if cancel_check is not None and cancel_check():
            raise PreprocessingCancelled("Route preparation cancelled")

    ledger_path = destination / "source_steps.sqlite.tmp"
    ledger_path.unlink(missing_ok=True)
    ledger = sqlite3.connect(ledger_path)
    counts: Counter[str] = Counter()
    distribution: Counter[str] = Counter()
    archive = destination / "routes.source.jsonl.gz"
    steps = destination / "route_steps.observations.jsonl.gz"
    try:
        ledger.execute("PRAGMA journal_mode=OFF")
        ledger.execute("PRAGMA synchronous=OFF")
        ledger.execute("""CREATE TABLE steps (
            source_id TEXT, reaction_hash TEXT, transformation_hash TEXT,
            reaction TEXT, line INTEGER, route_id TEXT, subtree_id TEXT,
            entry TEXT, occurrences INTEGER, PRIMARY KEY(source_id, reaction_hash))""")
        with (
            gzip.open(source, "rt", encoding="utf-8") as stream,
            _compressed_writer(archive) as out,
        ):
            for line_number, line in enumerate(stream, start=1):
                check_cancelled()
                route = json.loads(line)
                if not isinstance(route, dict) or not isinstance(
                    route.get("original_tree"), dict
                ):
                    raise ValueError(f"Invalid route source at line {line_number}")
                route_id = route.get("route_id")
                if not isinstance(route_id, str) or not route_id:
                    raise ValueError(f"Missing route ID at line {line_number}")
                expected = set(route["original_tree"]["reaction_ids"])
                seen: set[str] = set()
                out.write(
                    _json(
                        {
                            "schema_version": "released_route_source.v1",
                            "source": {
                                "source_file": source.name,
                                "source_file_sha256": source_hash,
                                "source_row_number": line_number,
                            },
                            "route_order": "unresolved",
                            "raw_route": route,
                        }
                    )
                    + "\n"
                )
                for subtree in route["subtrees"]:
                    subtree_id = subtree["subtree_id"]
                    prefix = f"{subtree_id}_"
                    for entry in subtree["reactions"]:
                        full_id = entry["_id"]
                        if not full_id.startswith(prefix):
                            raise ValueError(
                                f"Reaction ID prefix mismatch at line {line_number}"
                            )
                        source_id = full_id[len(prefix) :]
                        if source_id not in expected:
                            raise ValueError(
                                f"Unexpected original reaction ID at line {line_number}"
                            )
                        seen.add(source_id)
                        reaction = entry["reaction_smiles"]
                        if not isinstance(reaction, str):
                            raise ValueError(
                                f"Invalid reaction text at line {line_number}"
                            )
                        reaction_hash = _text_hash(reaction)
                        transformation, _ = _reaction_fields(reaction)
                        ledger.execute(
                            """INSERT INTO steps VALUES (?,?,?,?,?,?,?,?,1)
                            ON CONFLICT(source_id,reaction_hash) DO UPDATE
                            SET occurrences=occurrences+1""",
                            (
                                source_id,
                                reaction_hash,
                                _text_hash(transformation),
                                reaction,
                                line_number,
                                route_id,
                                subtree_id,
                                _json(entry),
                            ),
                        )
                        counts["subtree_reaction_occurrences"] += 1
                if seen != expected:
                    raise ValueError(
                        f"Missing original step structures at line {line_number}"
                    )
                counts["route_count"] += 1
                counts["original_step_occurrences"] += len(expected)
                distribution[str(route["original_tree"]["num_reactions"])] += 1
                if line_number % 10000 == 0:
                    ledger.commit()
                    emit(
                        "routes",
                        line_number,
                        f"Preserved {line_number:,} routes and their step memberships.",
                    )
        ledger.commit()
        counts["unique_source_reaction_ids"] = ledger.execute(
            "SELECT COUNT(DISTINCT source_id) FROM steps"
        ).fetchone()[0]
        ledger.execute(
            "CREATE TABLE conflicts AS SELECT source_id FROM steps "
            "GROUP BY source_id HAVING COUNT(*)>1"
        )
        ledger.execute("CREATE UNIQUE INDEX conflict_ids ON conflicts(source_id)")
        counts["conflicting_source_ids"] = ledger.execute(
            "SELECT COUNT(*) FROM conflicts"
        ).fetchone()[0]
        ledger.execute("CREATE TABLE csv_ids (source_id TEXT PRIMARY KEY)")
        with csv_source.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != ["id", "rxn_smiles"]:
                raise ValueError("Unexpected original reaction CSV schema")
            for row in reader:
                check_cancelled()
                identifier = row["id"]
                if not identifier.startswith("uspto_"):
                    raise ValueError("Unexpected original reaction CSV ID")
                source_id = identifier[len("uspto_") :]
                counts["original_csv_rows"] += 1
                try:
                    ledger.execute("INSERT INTO csv_ids VALUES (?)", (source_id,))
                except sqlite3.IntegrityError:
                    counts["duplicate_original_csv_ids"] += 1
                matched = ledger.execute(
                    "SELECT 1 FROM steps WHERE source_id=? "
                    "AND transformation_hash=? LIMIT 1",
                    (source_id, _text_hash(row["rxn_smiles"])),
                ).fetchone()
                if matched is None:
                    counts["original_csv_unmatched_rows"] += 1
                if counts["original_csv_rows"] % 100000 == 0:
                    ledger.commit()
                    emit(
                        "coverage",
                        counts["original_csv_rows"],
                        "Cross-checking original CSV IDs and structures "
                        "against route steps.",
                    )
        counts["route_ids_absent_from_csv"] = ledger.execute(
            "SELECT COUNT(DISTINCT source_id) FROM steps "
            "WHERE source_id NOT IN (SELECT source_id FROM csv_ids)"
        ).fetchone()[0]
        with _compressed_writer(steps) as out:
            cursor = ledger.execute("""SELECT s.source_id,s.reaction,s.line,
                s.route_id,s.subtree_id,
                s.entry,s.occurrences,c.source_id IS NOT NULL FROM steps s
                LEFT JOIN conflicts c ON c.source_id=s.source_id
                ORDER BY s.source_id,s.reaction_hash""")
            for (
                source_id,
                reaction,
                line,
                route_id,
                subtree_id,
                entry,
                occurrences,
                conflicting,
            ) in cursor:
                check_cancelled()
                observation = _step_observation(
                    source=source,
                    source_hash=source_hash,
                    source_id=source_id,
                    reaction=reaction,
                    line=line,
                    route_id=route_id,
                    subtree_id=subtree_id,
                    entry=json.loads(entry),
                    occurrences=occurrences,
                    conflicting=bool(conflicting),
                )
                out.write(_json(observation.to_dict()) + "\n")
                counts["output_step_observations"] += 1
                counts[f"step_status_{observation.ingestion_status}"] += 1
                if counts["output_step_observations"] % 10000 == 0:
                    emit(
                        "steps",
                        counts["output_step_observations"],
                        f"Exported {counts['output_step_observations']:,} "
                        "unique source-step observations.",
                    )
        coverage = not any(
            counts[key]
            for key in (
                "duplicate_original_csv_ids",
                "original_csv_unmatched_rows",
                "route_ids_absent_from_csv",
            )
        )
        report = {
            "definition_version": ROUTE_RELEASE_DEFINITION_VERSION,
            "source_path": str(source.resolve()),
            "source_sha256": source_hash,
            "original_csv_path": str(csv_source.resolve()),
            "original_csv_sha256": csv_hash,
            "counts": dict(sorted(counts.items())),
            "route_step_count_distribution": dict(sorted(distribution.items())),
            "coverage_complete": coverage,
            "reused": False,
            "output_sha256": {p.name: _sha256(p) for p in (archive, steps)},
            "limitations": [
                "Source-inferred routes; source array order is "
                "not experimental chronology.",
                "Supplied mapping is unvalidated; no reaction chemistry "
                "or registry resolution is performed.",
                "No train/test split or release chemistry-review status is assigned.",
                "CSV reaction export is superseded by route-derived steps "
                "retaining the middle agent field.",
            ],
        }
        temporary = manifest.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        temporary.replace(manifest)
        if not coverage:
            raise ValueError(f"Route/CSV coverage gate failed; inspect {manifest}")
        return report
    finally:
        ledger.close()
        ledger_path.unlink(missing_ok=True)
