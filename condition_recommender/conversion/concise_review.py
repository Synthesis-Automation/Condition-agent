"""Concise chemistry review export for canonical generic conversion records."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, Mapping, Optional

from .generic import GenericConversionCache, convert_record
from .input_schema import discover_csv_datasets, iter_csv_records

CONCISE_REACTION_REVIEW_SCHEMA_VERSION = "1.4"
CONCISE_REACTION_REVIEW_FIELDS = (
    "canonical_reaction_smiles",
    "reaction_display_label_detailed",
    "original_reaction_type",
    "detected_reaction_family",
    "detection_status",
    "transformation_class",
    "signature_id",
    "evidence_quality",
    "transformation_confidence",
    "reaction_completeness_status",
    "product_heavy_atom_coverage",
    "product_element_excess",
    "partial_transformation_class",
    "missing_product_atom_elements",
    "chemistry_status",
    "condition_status",
    "condition_stage_status",
    "index_eligibility",
    "admission_reasons",
    "warnings",
    "spectators",
    "steric_electronic_factors",
)

_SUBSCRIPT_TRANSLATION = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


@dataclass(frozen=True)
class ConciseReviewProgress:
    """Progress update for a recursive dataset-folder review export."""

    phase: str
    file_index: int
    file_count: int
    row_count: int
    current_file: str
    message: str


class ConciseReviewConversionCancelled(RuntimeError):
    """Raised when a caller cancels a folder conversion between records."""


def _readable_token(value: Any) -> str:
    return str(value or "").replace("_", " ").strip()


def _formula_text(value: Any) -> str:
    return str(value or "").translate(_SUBSCRIPT_TRANSLATION)


def _text_or_blank(value: Any) -> str:
    return "" if value is None else str(value)


def _enum_text(value: Any) -> str:
    return str(value.value if isinstance(value, Enum) else value or "")


def _spectator_summary(signature: Mapping[str, Any]) -> str:
    grouped: Dict[tuple[str, str], list[int]] = {}
    for group in signature.get("spectator_groups") or ():
        if not isinstance(group, Mapping):
            continue
        group_id = str(group.get("group_id") or "")
        label = str(group.get("chemist_label") or "")
        key = (label, group_id)
        distance = group.get("graph_distance")
        distances = grouped.setdefault(key, [])
        if distance is not None:
            distances.append(int(distance))
    values = []
    for (label, group_id), distances in sorted(grouped.items()):
        readable_group = _readable_token(group_id)
        display = _formula_text(label) or readable_group
        if readable_group and readable_group.casefold() != display.casefold():
            display += f" [{readable_group}]"
        count = max(1, len(distances))
        if count > 1:
            display = f"{count}× {display}"
        if distances:
            distance_text = "/".join(
                str(value) for value in sorted(set(distances))
            )
            display += f" (d={distance_text})"
        values.append(display)
    return "; ".join(values)


def _partner_label(partner: Mapping[str, Any]) -> str:
    role = _readable_token(partner.get("role"))
    if role:
        return role
    component = int(partner.get("component_index") or 0) + 1
    chemist_label = _formula_text(partner.get("chemist_label"))
    return f"P{component} ({chemist_label})" if chemist_label else f"P{component}"


def _steric_summary(steric: Mapping[str, Any]) -> str:
    steric_class = _readable_token(steric.get("class"))
    center_class = _readable_token(steric.get("center_substitution_class"))
    values = []
    if center_class and center_class == steric_class:
        values.append(f"{center_class} N center")
    else:
        if steric_class:
            values.append(steric_class)
        if center_class:
            values.append(f"{center_class} N center")
    ortho_count = steric.get("ortho_substituent_count")
    if ortho_count:
        values.append(f"ortho={int(ortho_count)}")
    attached = []
    for group in steric.get("attached_groups") or ():
        if not isinstance(group, Mapping):
            continue
        context = _readable_token(group.get("context")) or "group"
        attachment_class = _readable_token(
            group.get("attachment_carbon_class")
        )
        if attachment_class:
            text = f"{context} α-C {attachment_class}"
            if group.get("alpha_branched"):
                text += ", branched"
        else:
            text = context
        attached.append(text)
    values.extend(sorted(set(attached)))
    return ", ".join(values) or "unclassified"


def _electronic_summary(electronic: Mapping[str, Any]) -> str:
    electronic_class = _readable_token(electronic.get("class"))
    qualitative_sum = electronic.get("qualitative_sum")
    if qualitative_sum is None:
        return electronic_class or "unclassified"
    return f"{electronic_class or 'unclassified'} (q={float(qualitative_sum):+.2f})"


def _partner_environment_summary(signature: Mapping[str, Any]) -> str:
    values = []
    partners = sorted(
        (
            partner
            for partner in signature.get("partners") or ()
            if isinstance(partner, Mapping)
        ),
        key=lambda partner: (
            int(partner.get("component_index") or 0),
            str(partner.get("role") or ""),
            str(partner.get("partner_id") or ""),
        ),
    )
    for partner in partners:
        steric = partner.get("steric")
        electronic = partner.get("electronic")
        steric_value = steric if isinstance(steric, Mapping) else {}
        electronic_value = (
            electronic if isinstance(electronic, Mapping) else {}
        )
        if not steric_value and not electronic_value:
            continue
        values.append(
            f"{_partner_label(partner)}: "
            f"S={_steric_summary(steric_value)}; "
            f"E={_electronic_summary(electronic_value)}"
        )
    return " | ".join(values)


def iter_canonical_records(path: str | Path) -> Iterator[Dict[str, Any]]:
    """Stream records from canonical JSONL or a sharded manifest."""
    source = Path(path)
    if source.name.casefold() == "shard_manifest.json":
        payload = json.loads(source.read_text(encoding="utf-8"))
        if payload.get("artifact_type") != "generic_sharded_conversion":
            raise ValueError(f"Not a sharded conversion manifest: {source}")
        for entry in payload.get("shards") or ():
            if entry.get("status") != "complete":
                continue
            yield from iter_canonical_records(
                source.parent / str(entry["output_path"])
            )
        return
    handle = (
        gzip.open(source, mode="rt", encoding="utf-8")
        if source.suffix.casefold() == ".gz"
        else source.open(mode="r", encoding="utf-8")
    )
    with handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid record JSONL at {source}:{line_number}: {exc.msg}"
                ) from exc
            if not isinstance(value, dict):
                raise ValueError(
                    f"Record is not a JSON object at {source}:{line_number}"
                )
            yield value


def concise_reaction_review_row(record: Mapping[str, Any]) -> Dict[str, str]:
    """Select compact chemistry fields needed for rapid structural review."""
    display = record.get("reaction_display_label")
    display_value = display if isinstance(display, Mapping) else {}
    signature = record.get("reaction_signature")
    signature_value = signature if isinstance(signature, Mapping) else {}
    completeness = record.get("reaction_completeness")
    completeness_value = (
        completeness if isinstance(completeness, Mapping) else {}
    )
    partial = record.get("partial_product_transformation")
    partial_value = partial if isinstance(partial, Mapping) else {}
    warnings = sorted(
        {
            str(value)
            for value in (
                tuple(signature_value.get("warnings") or ())
                + tuple(completeness_value.get("warnings") or ())
                + tuple(display_value.get("warnings") or ())
                + tuple(partial_value.get("warnings") or ())
            )
            if value
        }
    )
    element_excess = completeness_value.get("product_element_excess")
    element_excess_value = (
        element_excess if isinstance(element_excess, Mapping) else {}
    )
    return {
        "canonical_reaction_smiles": str(
            record.get("canonical_reaction_smiles")
            or record.get("reaction_smiles")
            or ""
        ),
        "reaction_display_label_detailed": str(
            display_value.get("detailed") or ""
        ),
        "original_reaction_type": str(
            record.get("source_declared_family") or ""
        ),
        "detected_reaction_family": str(record.get("named_family") or ""),
        "detection_status": str(
            display_value.get("status")
            or record.get("reaction_label_status")
            or "unavailable"
        ),
        "transformation_class": str(
            record.get("transformation_class")
            or signature_value.get("transformation_class")
            or ""
        ),
        "signature_id": str(signature_value.get("signature_id") or ""),
        "evidence_quality": str(record.get("evidence_quality") or ""),
        "transformation_confidence": _text_or_blank(
            record.get("transformation_confidence")
        ),
        "reaction_completeness_status": str(
            completeness_value.get("status") or ""
        ),
        "product_heavy_atom_coverage": _text_or_blank(
            completeness_value.get("product_heavy_atom_coverage")
        ),
        "product_element_excess": (
            json.dumps(
                dict(sorted(element_excess_value.items())),
                separators=(",", ":"),
            )
            if element_excess_value
            else ""
        ),
        "partial_transformation_class": str(
            partial_value.get("transformation_class") or ""
        ),
        "missing_product_atom_elements": "; ".join(
            str(value)
            for value in partial_value.get("missing_product_atom_elements") or ()
        ),
        "chemistry_status": _enum_text(record.get("chemistry_status")),
        "condition_status": _enum_text(record.get("condition_status")),
        "condition_stage_status": _enum_text(
            record.get("condition_stage_status")
        ),
        "index_eligibility": _enum_text(record.get("index_eligibility")),
        "admission_reasons": "; ".join(
            str(value) for value in record.get("admission_reasons") or ()
        ),
        "warnings": "; ".join(warnings),
        "spectators": _spectator_summary(signature_value),
        "steric_electronic_factors": _partner_environment_summary(
            signature_value
        ),
    }


def export_concise_reaction_review_csv(
    records_path: str | Path,
    output_path: str | Path,
    *,
    progress_callback: Optional[Callable[[ConciseReviewProgress], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
    progress_interval: int = 1_000,
) -> Dict[str, Any]:
    """Write an Excel-friendly concise review CSV from canonical JSONL records."""
    if progress_interval < 1:
        raise ValueError("progress_interval must be positive")
    source = Path(records_path)
    if not source.is_file():
        raise ValueError(f"Canonical records file does not exist: {source}")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    row_count = 0
    completed = False
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=CONCISE_REACTION_REVIEW_FIELDS,
            )
            writer.writeheader()
            for record in iter_canonical_records(source):
                if cancel_check is not None and cancel_check():
                    raise ConciseReviewConversionCancelled(
                        "Review export cancelled"
                    )
                writer.writerow(concise_reaction_review_row(record))
                row_count += 1
                if (
                    progress_callback is not None
                    and row_count % progress_interval == 0
                ):
                    progress_callback(
                        ConciseReviewProgress(
                            phase="rows_exported",
                            file_index=1,
                            file_count=1,
                            row_count=row_count,
                            current_file=source.name,
                            message=(
                                f"Exported {row_count} review row(s)."
                            ),
                        )
                    )
        temporary.replace(destination)
        completed = True
    finally:
        if not completed and temporary.is_file():
            temporary.unlink()
    if progress_callback is not None:
        progress_callback(
            ConciseReviewProgress(
                phase="completed",
                file_index=1,
                file_count=1,
                row_count=row_count,
                current_file=source.name,
                message=f"Review CSV complete: {row_count} row(s).",
            )
        )
    return {
        "schema_version": CONCISE_REACTION_REVIEW_SCHEMA_VERSION,
        "artifact_type": "concise_reaction_review_csv",
        "records_path": str(source),
        "output_path": str(destination),
        "row_count": row_count,
        "columns": list(CONCISE_REACTION_REVIEW_FIELDS),
    }


def convert_dataset_folder_to_concise_review_csv(
    dataset_path: str | Path,
    output_path: str | Path,
    *,
    progress_callback: Optional[Callable[[ConciseReviewProgress], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
    progress_interval: int = 25,
) -> Dict[str, Any]:
    """Recursively convert source CSVs directly into the concise review view."""
    if progress_interval < 1:
        raise ValueError("progress_interval must be positive")
    source = Path(dataset_path)
    destination = Path(output_path)
    paths = tuple(
        path
        for path in discover_csv_datasets(source)
        if path.resolve() != destination.resolve()
    )
    if not paths:
        raise ValueError(f"No CSV datasets found at {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    row_count = 0
    cache = GenericConversionCache(max_entries=5_000)

    def report(
        phase: str,
        *,
        file_index: int,
        current_file: str,
        message: str,
    ) -> None:
        if progress_callback is not None:
            progress_callback(
                ConciseReviewProgress(
                    phase=phase,
                    file_index=file_index,
                    file_count=len(paths),
                    row_count=row_count,
                    current_file=current_file,
                    message=message,
                )
            )

    report(
        "discovered",
        file_index=0,
        current_file="",
        message=f"Found {len(paths)} CSV dataset file(s).",
    )
    completed = False
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=CONCISE_REACTION_REVIEW_FIELDS,
            )
            writer.writeheader()
            for file_index, path in enumerate(paths, start=1):
                if cancel_check is not None and cancel_check():
                    raise ConciseReviewConversionCancelled(
                        "Review conversion cancelled"
                    )
                relative_path = (
                    path.relative_to(source).as_posix()
                    if source.is_dir()
                    else path.name
                )
                report(
                    "file_started",
                    file_index=file_index,
                    current_file=relative_path,
                    message=(
                        f"Processing file {file_index}/{len(paths)}: "
                        f"{relative_path}"
                    ),
                )
                file_row_count = 0
                for raw_record in iter_csv_records(path):
                    if cancel_check is not None and cancel_check():
                        raise ConciseReviewConversionCancelled(
                            "Review conversion cancelled"
                        )
                    converted = convert_record(raw_record, cache=cache)
                    writer.writerow(
                        concise_reaction_review_row(converted.to_dict())
                    )
                    row_count += 1
                    file_row_count += 1
                    if file_row_count % progress_interval == 0:
                        report(
                            "rows_converted",
                            file_index=file_index,
                            current_file=relative_path,
                            message=(
                                f"{relative_path}: {file_row_count} row(s); "
                                f"{row_count} total."
                            ),
                        )
                report(
                    "file_completed",
                    file_index=file_index,
                    current_file=relative_path,
                    message=(
                        f"Completed {relative_path}: "
                        f"{file_row_count} reaction(s)."
                    ),
                )
        temporary.replace(destination)
        completed = True
    finally:
        if not completed and temporary.is_file():
            temporary.unlink()
    report(
        "completed",
        file_index=len(paths),
        current_file="",
        message=f"Finished {row_count} reaction(s).",
    )
    return {
        "schema_version": CONCISE_REACTION_REVIEW_SCHEMA_VERSION,
        "artifact_type": "concise_reaction_review_csv",
        "dataset_path": str(source),
        "output_path": str(destination),
        "source_file_count": len(paths),
        "row_count": row_count,
        "columns": list(CONCISE_REACTION_REVIEW_FIELDS),
    }


__all__ = [
    "CONCISE_REACTION_REVIEW_FIELDS",
    "CONCISE_REACTION_REVIEW_SCHEMA_VERSION",
    "ConciseReviewConversionCancelled",
    "ConciseReviewProgress",
    "concise_reaction_review_row",
    "convert_dataset_folder_to_concise_review_csv",
    "export_concise_reaction_review_csv",
    "iter_canonical_records",
]
