"""Streaming generic conversion engine for mixed reaction datasets."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import Any, Dict, Optional

from reactive_taxonomy import REACTION_SIGNATURE_SCHEMA_VERSION

from ..models import (
    AdmissionTier,
    ChemistryStatus,
    ConditionStatus,
    GENERIC_CONVERTER_DEFINITION_VERSION,
    IndexEligibility,
    OutcomeStatus,
)
from .flatten import GENERIC_REVIEW_FIELDS, flatten_generic_record
from .generic import convert_record
from .input_schema import discover_csv_datasets, iter_csv_records


def _markdown(report: Dict[str, Any]) -> str:
    return "\n".join(
        (
            "# Generic mixed-dataset conversion",
            "",
            f"- Source files: {report['file_count']}",
            f"- Source rows: {report['row_count']}",
            f"- Verified: {report['tier_counts']['verified']}",
            f"- Review: {report['tier_counts']['review']}",
            f"- Rejected: {report['tier_counts']['rejected']}",
            f"- Index eligible: {report['index_eligibility_counts']['eligible']}",
            f"- Signature coverage: {report['signature_count']}/{report['row_count']}",
            f"- Unique canonical reactions: {report['duplicate_summary']['unique_canonical_reactions']}",
            f"- Multi-recipe groups: {report['duplicate_summary']['multi_recipe_groups']}",
            "",
            "Nested `records.jsonl` is canonical. Tier CSV files are review views.",
            "",
        )
    )


def convert_datasets(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    max_rows: Optional[int] = None,
) -> Dict[str, Any]:
    """Convert one CSV or a directory of CSVs through the common engine."""
    paths = discover_csv_datasets(dataset_path)
    if not paths:
        raise ValueError(f"No CSV datasets found at {dataset_path}")
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    tier_counts = Counter()
    chemistry_quality_counts = Counter()
    condition_quality_counts = Counter()
    outcome_quality_counts = Counter()
    index_eligibility_counts = Counter()
    reason_counts = Counter()
    evidence_counts = Counter()
    transformation_counts = Counter()
    family_counts = Counter()
    reaction_scope_counts = Counter()
    source_counts = Counter()
    source_tiers: Dict[str, Counter[str]] = defaultdict(Counter)
    condition_status_counts = Counter()
    recipe_warning_counts = Counter()
    resolved_recipe_counts = Counter()
    role_bucket_counts = Counter()
    role_confidence_counts = Counter()
    canonical_groups: Counter[str] = Counter()
    group_recipes: Dict[str, set[str]] = defaultdict(set)
    row_count = signature_count = 0
    with ExitStack() as stack:
        jsonl = stack.enter_context(
            (destination / "records.jsonl").open("w", encoding="utf-8")
        )
        writers = {}
        for tier in AdmissionTier:
            handle = stack.enter_context(
                (destination / f"{tier.value}.csv").open(
                    "w", encoding="utf-8-sig", newline=""
                )
            )
            writer = csv.DictWriter(handle, fieldnames=GENERIC_REVIEW_FIELDS)
            writer.writeheader()
            writers[tier] = writer
        stop = False
        for path in paths:
            if stop:
                break
            for raw_record in iter_csv_records(path):
                if max_rows is not None and row_count >= max_rows:
                    stop = True
                    break
                record = convert_record(raw_record)
                payload = record.to_dict()
                jsonl.write(json.dumps(payload, ensure_ascii=False, sort_keys=True))
                jsonl.write("\n")
                writers[record.admission_tier].writerow(flatten_generic_record(record))
                row_count += 1
                tier = record.admission_tier.value
                tier_counts[tier] += 1
                chemistry_quality_counts[record.chemistry_status.value] += 1
                condition_quality_counts[record.condition_status.value] += 1
                outcome_quality_counts[record.outcome_status.value] += 1
                index_eligibility_counts[record.index_eligibility.value] += 1
                reason_counts.update(record.admission_reasons)
                evidence_counts[record.evidence_quality] += 1
                source_counts[record.source_dataset] += 1
                source_tiers[record.source_dataset][tier] += 1
                signature_count += int(record.reaction_signature is not None)
                if record.transformation_class:
                    transformation_counts[record.transformation_class] += 1
                if record.named_family:
                    family_counts[record.named_family] += 1
                topology = (record.reaction_signature or {}).get("topology") or {}
                reaction_scope = str(topology.get("reaction_scope") or "")
                if reaction_scope:
                    reaction_scope_counts[reaction_scope] += 1
                condition_status_counts.update(
                    record.condition_resolution.get("status_counts") or {}
                )
                recipe_warning_counts.update(
                    record.condition_resolution.get("recipe_warnings") or ()
                )
                if record.resolved_recipe_id:
                    resolved_recipe_counts[record.resolved_recipe_id] += 1
                for bucket in (
                    "catalysts",
                    "ligands",
                    "bases",
                    "acids",
                    "condensation_agents",
                    "oxidants",
                    "reductants",
                    "additives",
                    "solvents",
                    "other_components",
                ):
                    for component in (record.resolved_recipe or {}).get(bucket, ()):
                        role_bucket_counts[bucket] += 1
                        confidence = float(
                            component.get("primary_role_confidence", 0.0)
                        )
                        if confidence >= 0.9:
                            role_confidence_counts["high"] += 1
                        elif confidence >= 0.7:
                            role_confidence_counts["medium"] += 1
                        else:
                            role_confidence_counts["low"] += 1
                if record.canonical_reaction_id:
                    canonical_groups[record.canonical_reaction_id] += 1
                    group_recipes[record.canonical_reaction_id].add(
                        record.resolved_recipe_id or record.raw_recipe_id
                    )
    repeated_groups = sum(count > 1 for count in canonical_groups.values())
    multi_recipe_groups = sum(len(recipes) > 1 for recipes in group_recipes.values())
    report: Dict[str, Any] = {
        "schema_version": "1.2",
        "converter_version": GENERIC_CONVERTER_DEFINITION_VERSION,
        "reaction_signature_schema_version": REACTION_SIGNATURE_SCHEMA_VERSION,
        "dataset_path": str(Path(dataset_path)),
        "input_files": [str(path) for path in paths],
        "file_count": len(paths),
        "max_rows": max_rows,
        "row_count": row_count,
        "tier_counts": {tier.value: tier_counts[tier.value] for tier in AdmissionTier},
        "chemistry_status_counts": {
            status.value: chemistry_quality_counts[status.value]
            for status in ChemistryStatus
        },
        "condition_status_counts": {
            status.value: condition_quality_counts[status.value]
            for status in ConditionStatus
        },
        "outcome_status_counts": {
            status.value: outcome_quality_counts[status.value]
            for status in OutcomeStatus
        },
        "index_eligibility_counts": {
            status.value: index_eligibility_counts[status.value]
            for status in IndexEligibility
        },
        "reason_counts": dict(sorted(reason_counts.items())),
        "evidence_quality_counts": dict(sorted(evidence_counts.items())),
        "signature_count": signature_count,
        "signature_rate": round(signature_count / row_count, 6) if row_count else 0.0,
        "transformation_class_counts": dict(sorted(transformation_counts.items())),
        "named_family_counts": dict(sorted(family_counts.items())),
        "reaction_scope_counts": dict(sorted(reaction_scope_counts.items())),
        "condition_resolution_status_counts": dict(
            sorted(condition_status_counts.items())
        ),
        "resolved_recipe_count": len(resolved_recipe_counts),
        "recipe_warning_counts": dict(sorted(recipe_warning_counts.items())),
        "resolved_role_bucket_counts": dict(sorted(role_bucket_counts.items())),
        "role_confidence_counts": dict(sorted(role_confidence_counts.items())),
        "source_row_counts": dict(sorted(source_counts.items())),
        "source_tier_counts": {
            source: {tier.value: counts[tier.value] for tier in AdmissionTier}
            for source, counts in sorted(source_tiers.items())
        },
        "duplicate_summary": {
            "unique_canonical_reactions": len(canonical_groups),
            "repeated_reaction_groups": repeated_groups,
            "multi_recipe_groups": multi_recipe_groups,
        },
    }
    (destination / "conversion_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (destination / "conversion_report.md").write_text(
        _markdown(report), encoding="utf-8"
    )
    return report


__all__ = ["convert_datasets"]
