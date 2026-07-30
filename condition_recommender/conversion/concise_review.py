"""Concise chemistry review export for canonical generic conversion records."""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, Mapping, Optional

from reactive_taxonomy import render_reactivity_profile

from .generic import GenericConversionCache, convert_record
from .input_schema import discover_csv_datasets, iter_csv_records

CONCISE_REACTION_REVIEW_SCHEMA_VERSION = "2.3"
CONCISE_REACTION_REVIEW_FIELDS = (
    "canonical_reaction_smiles",
    "reaction_display_label_detailed",
    "original_reaction_type",
    "detected_reaction_family",
    "detection_status",
    "transformation_class",
    "signature_id",
    "fallback_descriptor_id",
    "fallback_evidence_mode",
    "fallback_retrieval_eligible",
    "fallback_ineligibility_reasons",
    "evidence_quality",
    "external_mapping_status",
    "external_mapping_provider",
    "external_mapping_confidence",
    "external_mapping_matched_hypotheses",
    "transformation_confidence",
    "stereochemical_changes",
    "reaction_completeness_status",
    "edit_hypothesis_count",
    "edit_hypothesis_ids",
    "product_heavy_atom_coverage",
    "product_element_excess",
    "partial_transformation_class",
    "removed_fragment",
    "installed_fragment",
    "installed_fragment_key",
    "fragment_source_status",
    "fragment_source_candidates",
    "missing_product_atom_elements",
    "chemistry_status",
    "condition_status",
    "condition_stage_status",
    "index_eligibility",
    "admission_reasons",
    "warnings",
    "spectators",
    "reactivity_profile",
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


def _stereochemical_change_summary(signature: Mapping[str, Any]) -> str:
    """Render compact, deterministic atom/bond stereo outcomes for CSV review."""
    summaries = []
    for change in signature.get("stereo_changes") or ():
        if not isinstance(change, Mapping):
            continue
        stereo_type = _readable_token(change.get("stereo_type")) or "stereo"
        old = str(change.get("old_descriptor") or "∅")
        new = str(change.get("new_descriptor") or "∅")
        outcome = _readable_token(change.get("change_type"))
        summary = f"{stereo_type}: {old}→{new}"
        if outcome:
            summary += f" ({outcome})"
        summaries.append(summary)
    return "; ".join(sorted(set(summaries)))


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
            distance_text = "/".join(str(value) for value in sorted(set(distances)))
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
        profile = partner.get("reactivity_profile")
        if not isinstance(profile, Mapping):
            continue
        values.append(
            f"{_partner_label(partner)}: {render_reactivity_profile(profile)}"
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
            yield from iter_canonical_records(source.parent / str(entry["output_path"]))
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
    fallback = record.get("fallback_descriptor")
    fallback_value = fallback if isinstance(fallback, Mapping) else {}
    completeness = record.get("reaction_completeness")
    completeness_value = completeness if isinstance(completeness, Mapping) else {}
    partial = record.get("partial_product_transformation")
    partial_value = partial if isinstance(partial, Mapping) else {}
    installed = partial_value.get("installed_fragment")
    installed_value = installed if isinstance(installed, Mapping) else {}
    hypotheses = tuple(
        value
        for value in record.get("reaction_edit_hypotheses") or ()
        if isinstance(value, Mapping)
    )
    external_mapping = record.get("external_atom_mapping")
    external_mapping_value = (
        external_mapping if isinstance(external_mapping, Mapping) else {}
    )
    warnings = sorted(
        {
            str(value)
            for value in (
                tuple(signature_value.get("warnings") or ())
                + tuple(completeness_value.get("warnings") or ())
                + tuple(display_value.get("warnings") or ())
                + tuple(partial_value.get("warnings") or ())
                + tuple(fallback_value.get("warnings") or ())
                + tuple(external_mapping_value.get("warnings") or ())
            )
            if value
        }
    )
    element_excess = completeness_value.get("product_element_excess")
    element_excess_value = element_excess if isinstance(element_excess, Mapping) else {}
    return {
        "canonical_reaction_smiles": str(
            record.get("canonical_reaction_smiles")
            or record.get("reaction_smiles")
            or ""
        ),
        "reaction_display_label_detailed": str(display_value.get("detailed") or ""),
        "original_reaction_type": str(record.get("source_declared_family") or ""),
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
        "fallback_descriptor_id": str(fallback_value.get("descriptor_id") or ""),
        "fallback_evidence_mode": str(fallback_value.get("evidence_mode") or ""),
        "fallback_retrieval_eligible": _text_or_blank(
            fallback_value.get("retrieval_eligible")
        ),
        "fallback_ineligibility_reasons": "; ".join(
            str(value) for value in fallback_value.get("ineligibility_reasons") or ()
        ),
        "evidence_quality": str(record.get("evidence_quality") or ""),
        "external_mapping_status": str(
            external_mapping_value.get("status") or ""
        ),
        "external_mapping_provider": str(
            (external_mapping_value.get("provider") or {}).get("provider_id")
            or ""
        ),
        "external_mapping_confidence": _text_or_blank(
            external_mapping_value.get("mapper_confidence")
        ),
        "external_mapping_matched_hypotheses": "; ".join(
            str(value)
            for value in external_mapping_value.get("matched_hypothesis_ids") or ()
        ),
        "transformation_confidence": _text_or_blank(
            record.get("transformation_confidence")
        ),
        "stereochemical_changes": _stereochemical_change_summary(signature_value),
        "reaction_completeness_status": str(completeness_value.get("status") or ""),
        "edit_hypothesis_count": str(len(hypotheses)),
        "edit_hypothesis_ids": "; ".join(
            str(value.get("hypothesis_id") or "")
            for value in hypotheses
            if value.get("hypothesis_id")
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
        "removed_fragment": str(partial_value.get("removed_fragment_smiles") or ""),
        "installed_fragment": str(
            installed_value.get("canonical_fragment_smiles") or ""
        ),
        "installed_fragment_key": str(installed_value.get("fragment_key") or ""),
        "fragment_source_status": str(installed_value.get("source_status") or ""),
        "fragment_source_candidates": json.dumps(
            installed_value.get("source_candidates") or (),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        if installed_value.get("source_candidates")
        else "",
        "missing_product_atom_elements": "; ".join(
            str(value)
            for value in partial_value.get("missing_product_atom_elements") or ()
        ),
        "chemistry_status": _enum_text(record.get("chemistry_status")),
        "condition_status": _enum_text(record.get("condition_status")),
        "condition_stage_status": _enum_text(record.get("condition_stage_status")),
        "index_eligibility": _enum_text(record.get("index_eligibility")),
        "admission_reasons": "; ".join(
            str(value) for value in record.get("admission_reasons") or ()
        ),
        "warnings": "; ".join(warnings),
        "spectators": _spectator_summary(signature_value),
        "reactivity_profile": _partner_environment_summary(signature_value),
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
                    raise ConciseReviewConversionCancelled("Review export cancelled")
                writer.writerow(concise_reaction_review_row(record))
                row_count += 1
                if progress_callback is not None and row_count % progress_interval == 0:
                    progress_callback(
                        ConciseReviewProgress(
                            phase="rows_exported",
                            file_index=1,
                            file_count=1,
                            row_count=row_count,
                            current_file=source.name,
                            message=(f"Exported {row_count} review row(s)."),
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
                        f"Processing file {file_index}/{len(paths)}: {relative_path}"
                    ),
                )
                file_row_count = 0
                for raw_record in iter_csv_records(path):
                    if cancel_check is not None and cancel_check():
                        raise ConciseReviewConversionCancelled(
                            "Review conversion cancelled"
                        )
                    converted = convert_record(raw_record, cache=cache)
                    writer.writerow(concise_reaction_review_row(converted.to_dict()))
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
                        f"Completed {relative_path}: {file_row_count} reaction(s)."
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
