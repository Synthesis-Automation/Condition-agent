"""Streaming quality and chemistry audit for mixed reaction datasets."""

from __future__ import annotations

import csv
import hashlib
import heapq
import json
from collections import Counter
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Tuple

from condition_registry.resolver import ConditionRegistry
from reactive_taxonomy import featurize_reaction

from .identities import canonical_reaction_identity, raw_recipe_id
from .input_schema import RawReactionRecord, discover_csv_datasets, iter_csv_records

_COVERAGE_FIELDS = (
    "reaction_id",
    "source_declared_family",
    "reaction_smiles",
    "yield_pct",
    "temperature_c",
    "time_h",
    "reference",
    "reactant_cas",
    "product_cas",
    "reagent_cas",
    "catalyst_cas",
    "solvent_cas",
    "experimental_procedure",
    "stages",
    "steps",
    "notes",
)


@dataclass
class _ReactionGroup:
    count: int = 0
    datasets: set[str] = field(default_factory=set)
    recipes: set[str] = field(default_factory=set)
    references: set[str] = field(default_factory=set)


def _present(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (tuple, list, dict, set)):
        return bool(value)
    return bool(str(value).strip())


def _schema_signature(path: Path) -> Tuple[str, ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
    return tuple(str(value).strip() for value in header)


def _sample_priority(record: RawReactionRecord) -> int:
    token = "\x1f".join(
        (record.source_dataset, record.reaction_id, record.reaction_smiles)
    )
    return int(hashlib.sha256(token.encode("utf-8")).hexdigest(), 16)


def _keep_sample(
    heap: List[Tuple[int, int, RawReactionRecord]],
    record: RawReactionRecord,
    sample_size: int,
    sequence: int,
) -> None:
    if sample_size <= 0:
        return
    priority = _sample_priority(record)
    item = (-priority, sequence, record)
    if len(heap) < sample_size:
        heapq.heappush(heap, item)
    elif priority < -heap[0][0]:
        heapq.heapreplace(heap, item)


def _counter_dict(counter: Counter[str]) -> Dict[str, int]:
    return dict(sorted(counter.items()))


def _audit_chemistry_sample(
    records: Iterable[RawReactionRecord],
) -> Dict[str, Any]:
    evidence = Counter()
    transformations = Counter()
    warnings = Counter()
    sampled = valid = signatures = named_families = 0
    for record in records:
        sampled += 1
        analysis = featurize_reaction(record.reaction_smiles)
        valid += int(analysis.valid)
        evidence[analysis.evidence_quality] += 1
        signatures += int(analysis.reaction_signature is not None)
        named_families += int(bool(analysis.named_family))
        if analysis.transformation_class:
            transformations[analysis.transformation_class] += 1
        warnings.update(analysis.warnings)
    return {
        "sampled_rows": sampled,
        "valid_reactions": valid,
        "signature_count": signatures,
        "named_family_count": named_families,
        "signature_rate": round(signatures / sampled, 6) if sampled else 0.0,
        "evidence_quality_counts": _counter_dict(evidence),
        "transformation_class_counts": _counter_dict(transformations),
        "warning_counts": _counter_dict(warnings),
    }


def _merge_chemistry_samples(reports: Iterable[Mapping[str, Any]]) -> Dict[str, Any]:
    evidence = Counter()
    transformations = Counter()
    warnings = Counter()
    sampled = valid = signatures = named_families = 0
    for report in reports:
        sampled += int(report["sampled_rows"])
        valid += int(report["valid_reactions"])
        signatures += int(report["signature_count"])
        named_families += int(report["named_family_count"])
        evidence.update(report["evidence_quality_counts"])
        transformations.update(report["transformation_class_counts"])
        warnings.update(report["warning_counts"])
    return {
        "sampled_rows": sampled,
        "valid_reactions": valid,
        "signature_count": signatures,
        "named_family_count": named_families,
        "signature_rate": round(signatures / sampled, 6) if sampled else 0.0,
        "evidence_quality_counts": _counter_dict(evidence),
        "transformation_class_counts": _counter_dict(transformations),
        "warning_counts": _counter_dict(warnings),
    }


def _registry_coverage(
    occurrences: Mapping[Tuple[str, str], int],
) -> Dict[str, Any]:
    registry = ConditionRegistry()
    by_field: Dict[str, Any] = {}
    for field_name in ("reagent_cas", "catalyst_cas", "solvent_cas"):
        identifiers = {
            identifier: count
            for (field, identifier), count in occurrences.items()
            if field == field_name
        }
        unique_statuses = Counter()
        occurrence_statuses = Counter()
        for identifier, count in identifiers.items():
            status = registry.resolve(cas=identifier).status
            unique_statuses[status] += 1
            occurrence_statuses[status] += count
        total = sum(occurrence_statuses.values())
        resolved = occurrence_statuses["resolved"]
        by_field[field_name] = {
            "unique_identifiers": len(identifiers),
            "unique_status_counts": _counter_dict(unique_statuses),
            "occurrences": total,
            "occurrence_status_counts": _counter_dict(occurrence_statuses),
            "occurrence_resolution_rate": round(resolved / total, 6)
            if total
            else 0.0,
        }
    return {
        "registry_size": registry.size,
        "by_source_field": by_field,
    }


def _duplicate_summary(groups: Mapping[str, _ReactionGroup]) -> Dict[str, Any]:
    repeated = multi_recipe = same_recipe = cross_dataset = 0
    pair_counts: Counter[str] = Counter()
    for group in groups.values():
        if group.count <= 1:
            continue
        repeated += 1
        if len(group.recipes) > 1:
            multi_recipe += 1
        else:
            same_recipe += 1
        if len(group.datasets) > 1:
            cross_dataset += 1
            for left, right in combinations(sorted(group.datasets), 2):
                pair_counts[f"{left} <> {right}"] += 1
    return {
        "unique_canonical_reactions": len(groups),
        "repeated_reaction_groups": repeated,
        "multi_recipe_reaction_groups": multi_recipe,
        "same_recipe_repeated_groups": same_recipe,
        "cross_dataset_reaction_groups": cross_dataset,
        "top_cross_dataset_overlaps": [
            {"datasets": pair, "reaction_groups": count}
            for pair, count in pair_counts.most_common(25)
        ],
    }


def _markdown(report: Mapping[str, Any]) -> str:
    duplicate = report["duplicate_groups"]
    chemistry = report["chemistry_sample"]
    lines = [
        "# Mixed reaction dataset audit",
        "",
        f"- Files: {report['file_count']}",
        f"- Rows: {report['row_count']}",
        f"- Canonicalizable reactions: {report['canonical_reaction_count']}",
        f"- Unique canonical reactions: {duplicate['unique_canonical_reactions']}",
        f"- Repeated reaction groups: {duplicate['repeated_reaction_groups']}",
        f"- Multi-recipe groups: {duplicate['multi_recipe_reaction_groups']}",
        f"- Cross-dataset groups: {duplicate['cross_dataset_reaction_groups']}",
        "",
        "## Chemistry sample",
        "",
        f"- Sampled rows: {chemistry['sampled_rows']}",
        f"- Valid reactions: {chemistry['valid_reactions']}",
        f"- Signatures: {chemistry['signature_count']}",
        f"- Signature rate: {chemistry['signature_rate']:.2%}",
        "",
        "## Field coverage",
        "",
        "| Field | Non-empty | Rate |",
        "|---|---:|---:|",
    ]
    for field_name, values in report["field_coverage"].items():
        lines.append(
            f"| `{field_name}` | {values['non_empty']} | {values['rate']:.2%} |"
        )
    lines.extend(
        [
            "",
            "Chemistry coverage is measured on the deterministic per-file sample; "
            "metadata, identity, registry, and duplicate metrics cover every row.",
            "",
        ]
    )
    return "\n".join(lines)


def audit_datasets(
    dataset_path: str | Path,
    output_dir: str | Path,
    *,
    chemistry_sample_per_file: int = 100,
) -> Dict[str, Any]:
    """Audit mixed CSV datasets while leaving all source files untouched."""
    paths = discover_csv_datasets(dataset_path)
    if not paths:
        raise ValueError(f"No CSV datasets found at {dataset_path}")
    field_counts = Counter()
    schema_counts = Counter()
    row_counts = Counter()
    yield_counts = Counter()
    condition_counts = Counter()
    identity_occurrences: Counter[Tuple[str, str]] = Counter()
    groups: Dict[str, _ReactionGroup] = {}
    chemistry_by_dataset: Dict[str, Dict[str, Any]] = {}
    total_rows = canonical_count = 0
    for path in paths:
        schema_counts["|".join(_schema_signature(path))] += 1
        sample_heap: List[Tuple[int, int, RawReactionRecord]] = []
        dataset_rows = 0
        for sequence, record in enumerate(iter_csv_records(path)):
            dataset_rows += 1
            total_rows += 1
            for field_name in _COVERAGE_FIELDS:
                if _present(getattr(record, field_name)):
                    field_counts[field_name] += 1
            if record.yield_pct is None:
                yield_counts["missing_or_non_numeric"] += 1
            elif 0.0 <= record.yield_pct <= 100.0:
                yield_counts["valid"] += 1
            else:
                yield_counts["out_of_range"] += 1
            condition_fields = (
                record.reagent_cas,
                record.catalyst_cas,
                record.solvent_cas,
            )
            condition_counts["any"] += int(any(condition_fields))
            condition_counts["all_three_source_fields"] += int(all(condition_fields))
            for field_name in ("reagent_cas", "catalyst_cas", "solvent_cas"):
                for identifier in getattr(record, field_name):
                    identity_occurrences[(field_name, identifier)] += 1
            identity = canonical_reaction_identity(record.reaction_smiles)
            if identity is not None:
                canonical_count += 1
                group = groups.setdefault(identity.reaction_id, _ReactionGroup())
                group.count += 1
                group.datasets.add(record.source_dataset)
                group.recipes.add(raw_recipe_id(record))
                if record.reference:
                    group.references.add(record.reference)
            _keep_sample(
                sample_heap,
                record,
                chemistry_sample_per_file,
                sequence,
            )
        row_counts[path.stem] = dataset_rows
        sampled_records = [item[2] for item in sample_heap]
        chemistry_by_dataset[path.stem] = _audit_chemistry_sample(sampled_records)
    field_coverage = {
        field_name: {
            "non_empty": field_counts[field_name],
            "rate": round(field_counts[field_name] / total_rows, 6),
        }
        for field_name in _COVERAGE_FIELDS
    }
    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "dataset_path": str(Path(dataset_path)),
        "file_count": len(paths),
        "row_count": total_rows,
        "row_counts_by_dataset": _counter_dict(row_counts),
        "schema_variants": [
            {"file_count": count, "columns": signature.split("|")}
            for signature, count in sorted(
                schema_counts.items(), key=lambda item: (-item[1], item[0])
            )
        ],
        "field_coverage": field_coverage,
        "yield_counts": _counter_dict(yield_counts),
        "condition_source_field_counts": _counter_dict(condition_counts),
        "canonical_reaction_count": canonical_count,
        "canonical_reaction_rate": round(canonical_count / total_rows, 6),
        "duplicate_groups": _duplicate_summary(groups),
        "condition_registry": _registry_coverage(identity_occurrences),
        "chemistry_sample_parameters": {
            "per_file": chemistry_sample_per_file,
            "selection": "lowest_sha256_source_identity",
        },
        "chemistry_sample": _merge_chemistry_samples(
            chemistry_by_dataset.values()
        ),
        "chemistry_sample_by_dataset": dict(sorted(chemistry_by_dataset.items())),
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "audit_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (destination / "audit_report.md").write_text(
        _markdown(report), encoding="utf-8"
    )
    return report


__all__ = ["audit_datasets"]
