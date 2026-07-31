"""Evaluate template-free reaction-core coverage on the sample reaction CSV.

Source reaction types are copied to the report for evaluation only.  They are
never passed to atom mapping, reaction featurization, center selection,
remote-subgraph classification, key construction, or generic-label rendering.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import replace
import json
from pathlib import Path
import statistics
import sys
import time
from typing import Any, Iterable, Mapping, Optional, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    EXTERNAL_MAPPING_EVIDENCE,
    ExternalAtomMappingResult,
    ReactionAnalysis,
    ReactionCoreProjection,
    RxnMapperProvider,
    featurize_reaction,
)


REPORT_SCHEMA_VERSION = "2.0"
DEFAULT_SOURCE = PROJECT_ROOT / "examples" / "sample_reactions.csv"
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "results" / "reaction_core_sample_report.csv"
)
DEFAULT_SUMMARY = (
    PROJECT_ROOT / "results" / "reaction_core_sample_report_summary.json"
)


REPORT_FIELDS = (
    "report_schema_version",
    "source_row_number",
    "reaction_type",
    "example",
    "rxn_smiles",
    "baseline_valid",
    "baseline_evidence_quality",
    "baseline_has_signature",
    "baseline_signature_id",
    "baseline_named_family",
    "baseline_transformation_class",
    "baseline_reaction_label",
    "baseline_reaction_label_status",
    "baseline_completeness_status",
    "baseline_warnings",
    "mapping_attempted",
    "mapping_status",
    "mapper_confidence",
    "reactant_mapping_coverage",
    "product_mapping_coverage",
    "shared_mapped_heavy_atom_count",
    "mapping_warnings",
    "mapping_error",
    "mapped_reaction_smiles",
    "mapped_analysis_has_signature",
    "core_available",
    "new_core_coverage",
    "core_id",
    "exact_core_key",
    "typed_core_key",
    "shape_core_key",
    "center_transition_key",
    "core_generic_label",
    "core_active_atom_count",
    "core_event_count",
    "core_primary_center_count",
    "core_atom_transitions",
    "remote_classes",
    "retained_remote_classes",
    "remote_fragments",
    "core_remote_subgraphs",
    "core_evidence",
    "core_confidence",
    "core_warnings",
    "review_priority",
    "review_reasons",
    "report_status",
)


def _read_rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return tuple(dict(row) for row in csv.DictReader(handle))


def _json(values: Any) -> str:
    return json.dumps(
        values,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _boolean(value: bool) -> str:
    return "true" if value else "false"


def _mapping_status(
    mapping: Optional[ExternalAtomMappingResult],
    *,
    baseline_valid: bool,
) -> str:
    if mapping is None:
        if not baseline_valid:
            return "not_attempted_invalid_reaction"
        return "not_attempted_supplied_mapping"
    if mapping.valid:
        return "valid"
    return "invalid"


def _mapped_analysis(
    mapping: ExternalAtomMappingResult,
    *,
    provider_id: str,
) -> Optional[ReactionAnalysis]:
    if (
        not mapping.valid
        or mapping.normalization is None
        or not mapping.mapped_reaction_smiles
    ):
        return None
    override = replace(
        mapping.normalization,
        evidence=EXTERNAL_MAPPING_EVIDENCE,
        confidence=mapping.normalization.confidence,
        warnings=tuple(
            sorted(
                set(mapping.normalization.warnings).union(
                    (
                        "EXTERNAL_MAPPING_WITHOUT_INTERNAL_CONSENSUS",
                        "EXTERNAL_MAPPING_REQUIRES_EXPERT_REVIEW",
                    )
                )
            )
        ),
    )
    return featurize_reaction(
        mapping.mapped_reaction_smiles,
        _mapped_edit_override=override,
        _mapped_provider=provider_id,
    )


def _core_atom_transitions(
    core: ReactionCoreProjection,
) -> list[dict[str, Any]]:
    return [
        {
            "atom_map_number": transition.atom_map_number,
            "before": (
                transition.before_state.concise_label
                if transition.before_state is not None
                else None
            ),
            "after": (
                transition.after_state.concise_label
                if transition.after_state is not None
                else None
            ),
            "incident_edit_count": transition.incident_edit_count,
            "stable_remote_subgraph_count": (
                transition.stable_remote_subgraph_count
            ),
            "role": transition.role,
        }
        for transition in core.atom_transitions
    ]


def _core_remote_subgraphs(
    core: ReactionCoreProjection,
) -> list[dict[str, Any]]:
    return [
        {
            "side": subgraph.side,
            "remote_class": subgraph.remote_class,
            "continuity": subgraph.continuity,
            "fragment_smiles": subgraph.fragment_smiles,
            "functional_group_ids": list(subgraph.functional_group_ids),
            "attachment_ports": [
                {
                    "core_atom_map_number": port.core_atom_map_number,
                    "attachment_atom_map_number": (
                        port.attachment_atom_map_number
                    ),
                    "attachment_element": port.attachment_element,
                    "bond_order": port.bond_order,
                }
                for port in subgraph.attachment_ports
            ],
        }
        for subgraph in core.remote_subgraphs
    ]


def _review(
    *,
    baseline: ReactionAnalysis,
    mapping: Optional[ExternalAtomMappingResult],
    mapped_analysis: Optional[ReactionAnalysis],
    core: Optional[ReactionCoreProjection],
) -> tuple[str, tuple[str, ...]]:
    reasons = []
    if not baseline.valid:
        reasons.append("invalid_reaction")
    if mapping is not None and not mapping.valid:
        reasons.append("external_mapping_failed")
    if core is None:
        reasons.append("reaction_core_unavailable")
    if mapping is not None and mapping.mapper_confidence is not None:
        if mapping.mapper_confidence < 0.5:
            reasons.append("low_mapper_confidence")
        elif mapping.mapper_confidence < 0.75:
            reasons.append("moderate_mapper_confidence")
    if core is not None:
        reasons.extend(core.warnings)
        if core.event_count > 1:
            reasons.append("multi_event_core")
        primary_count = sum(
            transition.role == "primary_center"
            for transition in core.atom_transitions
        )
        if primary_count > 2:
            reasons.append("many_selected_centers")
    if mapped_analysis is not None and mapped_analysis.reaction_signature is None:
        reasons.append("mapped_signature_unavailable")
    if baseline.reaction_signature is None:
        reasons.append("baseline_signature_unavailable")
    reasons = sorted(set(reasons))
    high_markers = {
        "external_mapping_failed",
        "invalid_reaction",
        "low_mapper_confidence",
        "reaction_core_unavailable",
        "REACTION_CORE_NO_OP_PRIMARY_CENTER",
        "REACTION_CORE_REMOTE_CONTINUITY_UNRESOLVED",
    }
    if high_markers.intersection(reasons):
        priority = "high"
    elif reasons:
        priority = "medium"
    else:
        priority = "low"
    return priority, tuple(reasons)


def _report_row(
    *,
    source_row_number: int,
    source: Mapping[str, str],
    baseline: ReactionAnalysis,
    mapping: Optional[ExternalAtomMappingResult],
    mapped_analysis: Optional[ReactionAnalysis],
) -> dict[str, Any]:
    core = (
        mapped_analysis.reaction_core
        if mapped_analysis is not None
        else baseline.reaction_core
    )
    completeness = baseline.reaction_completeness
    priority, review_reasons = _review(
        baseline=baseline,
        mapping=mapping,
        mapped_analysis=mapped_analysis,
        core=core,
    )
    mapping_attempted = mapping is not None
    mapping_warnings = mapping.warnings if mapping is not None else ()
    remote_classes = (
        sorted({subgraph.remote_class for subgraph in core.remote_subgraphs})
        if core is not None
        else []
    )
    retained_remote_classes = (
        sorted(
            {
                subgraph.remote_class
                for subgraph in core.remote_subgraphs
                if subgraph.continuity == "retained"
            }
        )
        if core is not None
        else []
    )
    remote_fragments = (
        sorted(
            {
                subgraph.fragment_smiles
                for subgraph in core.remote_subgraphs
                if subgraph.fragment_smiles
            }
        )
        if core is not None
        else []
    )
    baseline_signature = baseline.reaction_signature
    mapped_signature = (
        mapped_analysis.reaction_signature
        if mapped_analysis is not None
        else baseline_signature
    )
    report_status = (
        "core_available"
        if core is not None
        else "mapping_failed"
        if mapping is not None and not mapping.valid
        else "no_shared_mapped_center"
        if mapping is not None and mapping.valid
        else "core_unavailable"
    )
    return {
        "report_schema_version": REPORT_SCHEMA_VERSION,
        "source_row_number": source_row_number,
        "reaction_type": str(source.get("reaction_type") or ""),
        "example": str(source.get("example") or ""),
        "rxn_smiles": str(source.get("rxn_smiles") or ""),
        "baseline_valid": _boolean(baseline.valid),
        "baseline_evidence_quality": baseline.evidence_quality,
        "baseline_has_signature": _boolean(baseline_signature is not None),
        "baseline_signature_id": (
            baseline_signature.signature_id if baseline_signature else ""
        ),
        "baseline_named_family": baseline.named_family or "",
        "baseline_transformation_class": baseline.transformation_class or "",
        "baseline_reaction_label": baseline.reaction_label or "",
        "baseline_reaction_label_status": baseline.reaction_label_status,
        "baseline_completeness_status": (
            completeness.status if completeness is not None else ""
        ),
        "baseline_warnings": _json(list(baseline.warnings)),
        "mapping_attempted": _boolean(mapping_attempted),
        "mapping_status": _mapping_status(
            mapping,
            baseline_valid=baseline.valid,
        ),
        "mapper_confidence": (
            round(float(mapping.mapper_confidence), 6)
            if mapping is not None and mapping.mapper_confidence is not None
            else ""
        ),
        "reactant_mapping_coverage": (
            round(mapping.reactant_mapping_coverage, 6)
            if mapping is not None
            else ""
        ),
        "product_mapping_coverage": (
            round(mapping.product_mapping_coverage, 6)
            if mapping is not None
            else ""
        ),
        "shared_mapped_heavy_atom_count": (
            mapping.shared_mapped_heavy_atom_count
            if mapping is not None
            else ""
        ),
        "mapping_warnings": _json(list(mapping_warnings)),
        "mapping_error": (
            str(mapping.error or "") if mapping is not None else ""
        ),
        "mapped_reaction_smiles": (
            str(mapping.mapped_reaction_smiles or "")
            if mapping is not None
            else ""
        ),
        "mapped_analysis_has_signature": _boolean(mapped_signature is not None),
        "core_available": _boolean(core is not None),
        "new_core_coverage": _boolean(
            core is not None and baseline_signature is None
        ),
        "core_id": core.core_id if core is not None else "",
        "exact_core_key": core.exact_core_key if core is not None else "",
        "typed_core_key": core.typed_core_key if core is not None else "",
        "shape_core_key": core.shape_core_key if core is not None else "",
        "center_transition_key": (
            core.center_transition_key if core is not None else ""
        ),
        "core_generic_label": core.generic_label if core is not None else "",
        "core_active_atom_count": (
            core.active_atom_count if core is not None else ""
        ),
        "core_event_count": core.event_count if core is not None else "",
        "core_primary_center_count": (
            sum(
                transition.role == "primary_center"
                for transition in core.atom_transitions
            )
            if core is not None
            else ""
        ),
        "core_atom_transitions": (
            _json(_core_atom_transitions(core)) if core is not None else "[]"
        ),
        "remote_classes": _json(remote_classes),
        "retained_remote_classes": _json(retained_remote_classes),
        "remote_fragments": _json(remote_fragments),
        "core_remote_subgraphs": (
            _json(_core_remote_subgraphs(core)) if core is not None else "[]"
        ),
        "core_evidence": core.evidence if core is not None else "",
        "core_confidence": (
            round(core.confidence, 6) if core is not None else ""
        ),
        "core_warnings": (
            _json(list(core.warnings)) if core is not None else "[]"
        ),
        "review_priority": priority,
        "review_reasons": _json(list(review_reasons)),
        "report_status": report_status,
    }


def _counts(values: Iterable[str]) -> dict[str, int]:
    return dict(sorted(Counter(values).items()))


def _summary(
    report_rows: Sequence[Mapping[str, Any]],
    *,
    source: Path,
    elapsed_seconds: float,
) -> dict[str, Any]:
    total = len(report_rows)
    core_rows = [
        row for row in report_rows if row["core_available"] == "true"
    ]
    mapped_rows = [
        row for row in report_rows if row["mapping_attempted"] == "true"
    ]
    valid_mapping_rows = [
        row for row in mapped_rows if row["mapping_status"] == "valid"
    ]
    confidence_values = [
        float(row["mapper_confidence"])
        for row in valid_mapping_rows
        if row["mapper_confidence"] != ""
    ]
    clusters: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    shape_clusters: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in core_rows:
        clusters[str(row["center_transition_key"])].append(row)
        shape_clusters[str(row["shape_core_key"])].append(row)
    repeated_clusters = {
        key: rows for key, rows in clusters.items() if len(rows) > 1
    }
    mixed_clusters = {
        key: rows
        for key, rows in repeated_clusters.items()
        if len({str(row["reaction_type"]) for row in rows}) > 1
    }
    repeated_shape_clusters = {
        key: rows for key, rows in shape_clusters.items() if len(rows) > 1
    }
    mixed_shape_clusters = {
        key: rows
        for key, rows in repeated_shape_clusters.items()
        if len({str(row["reaction_type"]) for row in rows}) > 1
    }
    by_type = {}
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in report_rows:
        grouped[str(row["reaction_type"])].append(row)
    for reaction_type, rows in sorted(grouped.items()):
        by_type[reaction_type] = {
            "rows": len(rows),
            "baseline_signatures": sum(
                row["baseline_has_signature"] == "true" for row in rows
            ),
            "valid_mappings": sum(
                row["mapping_status"] == "valid" for row in rows
            ),
            "cores": sum(row["core_available"] == "true" for row in rows),
            "new_core_coverage": sum(
                row["new_core_coverage"] == "true" for row in rows
            ),
            "high_review_priority": sum(
                row["review_priority"] == "high" for row in rows
            ),
        }
    return {
        "report_schema_version": REPORT_SCHEMA_VERSION,
        "source": str(source),
        "rows": total,
        "elapsed_seconds": round(elapsed_seconds, 3),
        "baseline": {
            "valid": sum(row["baseline_valid"] == "true" for row in report_rows),
            "with_signature": sum(
                row["baseline_has_signature"] == "true"
                for row in report_rows
            ),
            "with_named_family": sum(
                bool(row["baseline_named_family"]) for row in report_rows
            ),
            "evidence_quality": _counts(
                str(row["baseline_evidence_quality"])
                for row in report_rows
            ),
            "completeness_status": _counts(
                str(row["baseline_completeness_status"])
                for row in report_rows
            ),
        },
        "mapping": {
            "attempted": len(mapped_rows),
            "valid": len(valid_mapping_rows),
            "invalid": len(mapped_rows) - len(valid_mapping_rows),
            "median_confidence": (
                round(statistics.median(confidence_values), 6)
                if confidence_values
                else None
            ),
            "mean_confidence": (
                round(statistics.mean(confidence_values), 6)
                if confidence_values
                else None
            ),
        },
        "reaction_core": {
            "available": len(core_rows),
            "unavailable": total - len(core_rows),
            "new_coverage_without_baseline_signature": sum(
                row["new_core_coverage"] == "true"
                for row in report_rows
            ),
            "unique_center_transition_keys": len(clusters),
            "repeated_center_transition_keys": len(repeated_clusters),
            "rows_in_repeated_center_clusters": sum(
                len(rows) for rows in repeated_clusters.values()
            ),
            "mixed_source_label_clusters": len(mixed_clusters),
            "rows_in_mixed_source_label_clusters": sum(
                len(rows) for rows in mixed_clusters.values()
            ),
            "unique_shape_core_keys": len(shape_clusters),
            "repeated_shape_core_keys": len(repeated_shape_clusters),
            "rows_in_repeated_shape_clusters": sum(
                len(rows) for rows in repeated_shape_clusters.values()
            ),
            "mixed_source_label_shape_clusters": len(mixed_shape_clusters),
            "rows_in_mixed_source_label_shape_clusters": sum(
                len(rows) for rows in mixed_shape_clusters.values()
            ),
            "event_count": _counts(
                str(row["core_event_count"]) for row in core_rows
            ),
            "primary_center_count": _counts(
                str(row["core_primary_center_count"]) for row in core_rows
            ),
            "generic_labels": _counts(
                str(row["core_generic_label"]) for row in core_rows
            ),
            "remote_classes": _counts(
                remote_class
                for row in core_rows
                for remote_class in json.loads(
                    str(row["remote_classes"])
                )
            ),
        },
        "review_priority": _counts(
            str(row["review_priority"]) for row in report_rows
        ),
        "review_reasons": _counts(
            reason
            for row in report_rows
            for reason in json.loads(str(row["review_reasons"]))
        ),
        "report_status": _counts(
            str(row["report_status"]) for row in report_rows
        ),
        "by_reaction_type": by_type,
    }


def build_report(
    source: Path,
    output: Path,
    summary_output: Path,
    *,
    batch_size: int,
    outer_batch_size: int,
) -> dict[str, Any]:
    """Build and persist the sample reaction-core evaluation report."""
    started = time.perf_counter()
    source_rows = _read_rows(source)
    reactions = tuple(str(row.get("rxn_smiles") or "").strip() for row in source_rows)
    unique_reactions = tuple(dict.fromkeys(reaction for reaction in reactions if reaction))
    print(
        f"Loaded {len(source_rows)} rows ({len(unique_reactions)} unique reactions)",
        flush=True,
    )
    baseline_by_reaction = {
        reaction: featurize_reaction(reaction) for reaction in unique_reactions
    }
    to_map = tuple(
        reaction
        for reaction in unique_reactions
        if baseline_by_reaction[reaction].valid
        and not any(
            component.atom_mapped
            for component in (
                *baseline_by_reaction[reaction].reactants,
                *baseline_by_reaction[reaction].products,
            )
        )
    )
    print(
        f"Baseline complete; mapping {len(to_map)} unique unmapped reactions",
        flush=True,
    )
    if to_map and not RxnMapperProvider.is_available():
        raise RuntimeError(
            "RXNMapper is unavailable; install requirements-mapping.txt"
        )
    provider = RxnMapperProvider(batch_size=batch_size)
    mappings: dict[str, ExternalAtomMappingResult] = {}
    for offset in range(0, len(to_map), outer_batch_size):
        batch = to_map[offset : offset + outer_batch_size]
        results = provider.map_reactions(batch)
        mappings.update(zip(batch, results))
        print(
            f"Mapped {min(offset + len(batch), len(to_map))}/{len(to_map)}",
            flush=True,
        )
    mapped_analyses: dict[str, Optional[ReactionAnalysis]] = {}
    for index, reaction in enumerate(to_map, start=1):
        mapped_analyses[reaction] = _mapped_analysis(
            mappings[reaction],
            provider_id=provider.metadata.provider_id,
        )
        if index % 100 == 0 or index == len(to_map):
            print(
                f"Built mapped analyses {index}/{len(to_map)}",
                flush=True,
            )
    report_rows = [
        _report_row(
            source_row_number=index,
            source=row,
            baseline=baseline_by_reaction[reaction],
            mapping=mappings.get(reaction),
            mapped_analysis=mapped_analyses.get(reaction),
        )
        for index, (row, reaction) in enumerate(
            zip(source_rows, reactions),
            start=2,
        )
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS)
        writer.writeheader()
        writer.writerows(report_rows)
    summary = _summary(
        report_rows,
        source=source,
        elapsed_seconds=time.perf_counter() - started,
    )
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(
        json.dumps(summary, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {output}", flush=True)
    print(f"Wrote {summary_output}", flush=True)
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--outer-batch-size", type=int, default=64)
    args = parser.parse_args(argv)
    summary = build_report(
        args.source,
        args.output,
        args.summary,
        batch_size=args.batch_size,
        outer_batch_size=args.outer_batch_size,
    )
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
