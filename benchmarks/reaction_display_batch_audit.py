"""Render display-minimized graphs for every row in a reaction review CSV."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any, Dict, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    build_reaction_display_projection,
    featurize_reaction,
    reaction_render_context_from_analysis,
)
from visualization import build_reaction_display_graphic  # noqa: E402


DEFAULT_SOURCE = PROJECT_ROOT / "datasets" / "literature" / "reaction_review.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "reaction_review_core_graphs"
FAILURE_COLUMNS = (
    "source_data_row",
    "source_csv_line",
    "canonical_reaction_smiles",
    "reaction_label",
    "original_reaction_type",
    "source_reaction_core_status",
    "source_reaction_core_quality_status",
    "failure_stage",
    "exception_type",
    "failure_reason",
    "analysis_valid",
    "computed_evidence_status",
    "computed_warnings",
)


@dataclass(frozen=True)
class _AuditOutcome:
    """Cached graph-generation result for one unique reaction SMILES."""

    success: bool
    image_bytes: bytes | None
    failure_stage: str
    exception_type: str
    failure_reason: str
    analysis_valid: bool | None
    evidence_status: str
    warnings: tuple[str, ...]
    definition_id: str


def _failure(
    *,
    stage: str,
    reason: str,
    exception: Exception | None = None,
    analysis_valid: bool | None = None,
    evidence_status: str = "",
    warnings: tuple[str, ...] = (),
) -> _AuditOutcome:
    return _AuditOutcome(
        success=False,
        image_bytes=None,
        failure_stage=stage,
        exception_type=type(exception).__name__ if exception else "",
        failure_reason=reason,
        analysis_valid=analysis_valid,
        evidence_status=evidence_status,
        warnings=warnings,
        definition_id="",
    )


def _audit_reaction(reaction_smiles: str) -> _AuditOutcome:
    if not reaction_smiles.strip():
        return _failure(stage="input", reason="missing reaction SMILES")
    try:
        analysis = featurize_reaction(reaction_smiles)
    except Exception as exc:
        return _failure(
            stage="featurization",
            reason=str(exc),
            exception=exc,
        )
    analysis_valid = bool(getattr(analysis, "valid", False))
    analysis_warnings = tuple(str(value) for value in analysis.warnings)
    if not analysis_valid:
        return _failure(
            stage="analysis_validation",
            reason="reaction analysis is invalid",
            analysis_valid=False,
            warnings=analysis_warnings,
        )
    core = analysis.reaction_core
    if core is None:
        return _failure(
            stage="reaction_core",
            reason="reaction core is unavailable",
            analysis_valid=True,
            warnings=analysis_warnings,
        )
    evidence_status = str(core.evidence_status)
    try:
        projection = build_reaction_display_projection(
            reaction_render_context_from_analysis(analysis)
        )
    except Exception as exc:
        return _failure(
            stage="display_projection",
            reason=str(exc),
            exception=exc,
            analysis_valid=True,
            evidence_status=evidence_status,
            warnings=analysis_warnings + tuple(core.warnings),
        )
    try:
        graphic = build_reaction_display_graphic(
            projection,
            size=(1100, 280),
            image_format="svg",
        )
    except Exception as exc:
        return _failure(
            stage="graphic_rendering",
            reason=str(exc),
            exception=exc,
            analysis_valid=True,
            evidence_status=evidence_status,
            warnings=projection.warnings,
        )
    return _AuditOutcome(
        success=True,
        image_bytes=graphic.image_bytes,
        failure_stage="",
        exception_type="",
        failure_reason="",
        analysis_valid=True,
        evidence_status=evidence_status,
        warnings=projection.warnings,
        definition_id=projection.definition_id,
    )


def _graph_filename(data_row: int, reaction_smiles: str) -> str:
    digest = hashlib.sha256(reaction_smiles.encode("utf-8")).hexdigest()[:12]
    return f"row_{data_row:06d}_{digest}.svg"


def _relative_or_absolute(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def build_batch_audit(
    *,
    source_path: Path,
    output_dir: Path,
) -> Dict[str, Any]:
    """Render every CSV row and write structured failure and summary reports."""
    started = time.perf_counter()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"output directory is not empty: {output_dir.resolve()}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    graph_dir = output_dir / "graphs"
    graph_dir.mkdir()

    with source_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("reaction review CSV has no data rows")
    if "canonical_reaction_smiles" not in rows[0]:
        raise ValueError(
            "reaction review CSV requires canonical_reaction_smiles"
        )

    cache: Dict[str, _AuditOutcome] = {}
    unsuccessful_rows = []
    successful_rows = 0
    definition_ids = set()
    failure_by_stage: Dict[str, int] = {}
    for data_row, row in enumerate(rows, start=1):
        reaction_smiles = str(row.get("canonical_reaction_smiles") or "").strip()
        outcome = cache.get(reaction_smiles)
        if outcome is None:
            outcome = _audit_reaction(reaction_smiles)
            cache[reaction_smiles] = outcome
        if outcome.success:
            if outcome.image_bytes is None:
                raise AssertionError("successful audit outcome has no graphic")
            filename = _graph_filename(data_row, reaction_smiles)
            (graph_dir / filename).write_bytes(outcome.image_bytes)
            successful_rows += 1
            if outcome.definition_id:
                definition_ids.add(outcome.definition_id)
            continue
        failure_by_stage[outcome.failure_stage] = (
            failure_by_stage.get(outcome.failure_stage, 0) + 1
        )
        unsuccessful_rows.append(
            {
                "source_data_row": data_row,
                "source_csv_line": data_row + 1,
                "canonical_reaction_smiles": reaction_smiles,
                "reaction_label": str(row.get("reaction_label") or ""),
                "original_reaction_type": str(
                    row.get("original_reaction_type") or ""
                ),
                "source_reaction_core_status": str(
                    row.get("reaction_core_status") or ""
                ),
                "source_reaction_core_quality_status": str(
                    row.get("reaction_core_quality_status") or ""
                ),
                "failure_stage": outcome.failure_stage,
                "exception_type": outcome.exception_type,
                "failure_reason": outcome.failure_reason,
                "analysis_valid": (
                    "" if outcome.analysis_valid is None else outcome.analysis_valid
                ),
                "computed_evidence_status": outcome.evidence_status,
                "computed_warnings": "; ".join(sorted(set(outcome.warnings))),
            }
        )

    failure_path = output_dir / "unsuccessful_rows.csv"
    with failure_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FAILURE_COLUMNS)
        writer.writeheader()
        writer.writerows(unsuccessful_rows)

    unique_successes = sum(value.success for value in cache.values())
    summary: Dict[str, Any] = {
        "artifact_type": "reaction_display_batch_audit",
        "source_csv": _relative_or_absolute(source_path),
        "output_directory": _relative_or_absolute(output_dir),
        "graph_directory": _relative_or_absolute(graph_dir),
        "failure_csv": _relative_or_absolute(failure_path),
        "input_rows": len(rows),
        "unique_reactions": len(cache),
        "successful_rows": successful_rows,
        "unsuccessful_rows": len(unsuccessful_rows),
        "success_rate": successful_rows / len(rows),
        "unique_successful_reactions": unique_successes,
        "unique_unsuccessful_reactions": len(cache) - unique_successes,
        "failure_by_stage": dict(sorted(failure_by_stage.items())),
        "projection_definition_ids": sorted(definition_ids),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }
    (output_dir / "batch_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    summary = build_batch_audit(
        source_path=args.source,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
