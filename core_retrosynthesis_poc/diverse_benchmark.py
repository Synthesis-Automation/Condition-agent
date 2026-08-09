"""Leakage-safe benchmark across structurally diverse reaction archetypes."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Sequence

from retrosynthesis_poc.library import iter_rows

from .comparison import split_by_reference
from .ensemble import merge_baseline_first
from .generic_compiler import compile_generic_templates
from .generic_library import (
    build_generic_library,
    save_generic_library,
)
from .generic_search import disconnect_generic_target


DIVERSE_COHORTS = {
    "acyl_substitution": ("Amide_formation*.jsonl.gz",),
    "c_c_coupling": ("suzuki_miyaura*.jsonl.gz",),
    "carbonyl_oxidation": ("DessMartin_periodinane*.jsonl.gz",),
    "carbonyl_reduction": ("NaBH4_carbonyl_reductions*.jsonl.gz",),
    "conjugate_addition": ("Michael_addition*.jsonl.gz",),
    "carbonyl_condensation": ("Reductive_amination*.jsonl.gz",),
    "ring_formation": ("CuAAC_azidealkyne*.jsonl.gz",),
}


def load_diverse_rows(
    source: str | Path,
    *,
    max_rows_per_cohort: int = 200,
) -> tuple[Dict[str, Any], ...]:
    """Sample named cohorts while retaining structural routing independence."""

    if max_rows_per_cohort < 1:
        raise ValueError("max rows per cohort must be positive")
    values = []
    for cohort, patterns in DIVERSE_COHORTS.items():
        count = 0
        for row in iter_rows(source, include=patterns):
            if count >= max_rows_per_cohort:
                break
            values.append({**row, "_sampling_cohort": cohort})
            count += 1
    return tuple(values)


def _empty_metrics() -> Dict[str, Any]:
    return {
        "targets": 0,
        "targets_with_candidates": 0,
        "candidate_count": 0,
        "valid_candidate_count": 0,
        "top1_exact_precursor": 0,
        "top5_exact_precursor": 0,
        "top10_exact_precursor": 0,
        "top1_center": 0,
        "top5_center": 0,
        "top10_center": 0,
    }


def _record(
    metrics: Dict[str, Any],
    candidates: Sequence[Any],
    expected_precursors: str,
) -> tuple[int | None, int | None]:
    metrics["targets"] += 1
    metrics["targets_with_candidates"] += int(bool(candidates))
    metrics["candidate_count"] += len(candidates)
    metrics["valid_candidate_count"] += sum(
        candidate.forward_validation_status not in {"invalid", "unresolved"}
        for candidate in candidates
    )
    values = [candidate.precursor_smiles for candidate in candidates]
    exact_rank = values.index(expected_precursors) + 1 if expected_precursors in values else None
    archetype_rank = 1 if candidates else None
    for limit in (1, 5, 10):
        metrics[f"top{limit}_exact_precursor"] += int(
            exact_rank is not None and exact_rank <= limit
        )
        metrics[f"top{limit}_center"] += int(bool(candidates[:limit]))
    return exact_rank, archetype_rank


def _finalize(metrics: Dict[str, Any]) -> Dict[str, Any]:
    targets = max(1, int(metrics["targets"]))
    candidates = max(1, int(metrics["candidate_count"]))
    value = dict(metrics)
    value["target_coverage"] = metrics["targets_with_candidates"] / targets
    value["mean_candidates_per_target"] = metrics["candidate_count"] / targets
    value["valid_candidate_fraction"] = metrics["valid_candidate_count"] / candidates
    for limit in (1, 5, 10):
        value[f"top{limit}_exact_precursor_recall"] = (
            metrics[f"top{limit}_exact_precursor"] / targets
        )
        value[f"top{limit}_center_recall"] = (
            metrics[f"top{limit}_center"] / targets
        )
    return value


def _balanced_targets(
    cases: Sequence[tuple[Dict[str, Any], str, str, str]],
    max_targets_per_transformation: int,
) -> tuple[tuple[Dict[str, Any], str, str, str], ...]:
    grouped: dict[str, list[tuple[Dict[str, Any], str, str, str]]] = {}
    for case in cases:
        grouped.setdefault(case[3], []).append(case)
    selected = []
    for transformation, values in sorted(grouped.items()):
        ranked = sorted(
            values,
            key=lambda case: hashlib.sha256(
                f"{transformation}\0{case[0].get('reaction_id')}".encode()
            ).digest(),
        )
        selected.extend(ranked[:max_targets_per_transformation])
    return tuple(selected)


def run_diverse_benchmark(
    rows: Sequence[Dict[str, Any]],
    output_directory: str | Path,
    *,
    test_fraction: float = 0.2,
    max_targets_per_transformation: int = 10,
    top_k: int = 10,
    max_candidates_to_validate: int = 30,
) -> Dict[str, Any]:
    """Build baseline/core libraries and evaluate held-out diverse targets."""

    train, test = split_by_reference(rows, test_fraction=test_fraction)
    baseline_library = build_generic_library(train, engine="rdchiral")
    core_library = build_generic_library(train, engine="reaction_core")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    save_generic_library(baseline_library, output / "baseline_templates.json.gz")
    save_generic_library(core_library, output / "core_templates.json.gz")

    cases = []
    eligibility = Counter()
    for row in test:
        core = compile_generic_templates(row, engine="reaction_core", levels=("L1",))
        baseline = compile_generic_templates(row, engine="rdchiral")
        template = core.templates[0] if core.templates else (
            baseline.templates[0] if baseline.templates else None
        )
        if template is None:
            eligibility[str(core.rejection_reason or baseline.rejection_reason)] += 1
            continue
        precedent = template.precedents[0]
        cases.append(
            (
                row,
                precedent.product_smiles,
                precedent.precursor_smiles,
                template.transformation_kind,
            )
        )
    selected_cases = _balanced_targets(cases, max_targets_per_transformation)
    methods = (
        "baseline",
        "core_l1_context",
        "core_l2_context",
        "core_context",
        "ensemble_baseline_core_context",
    )
    metrics = {method: _empty_metrics() for method in methods}
    metrics_by_kind: dict[str, dict[str, Dict[str, Any]]] = {}
    target_results = []
    for row, target, expected, transformation in selected_cases:
        baseline = disconnect_generic_target(
            target,
            baseline_library,
            transformations=(transformation,),
            top_k=top_k,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        l1 = disconnect_generic_target(
            target,
            core_library,
            transformations=(transformation,),
            levels=("L1",),
            top_k=top_k,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        l2 = disconnect_generic_target(
            target,
            core_library,
            transformations=(transformation,),
            levels=("L2",),
            top_k=top_k,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        core_combined = disconnect_generic_target(
            target,
            core_library,
            transformations=(transformation,),
            levels=("L1", "L2"),
            top_k=top_k,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        ensemble = merge_baseline_first(baseline, core_combined, top_k)
        candidate_sets = {
            "baseline": baseline,
            "core_l1_context": l1,
            "core_l2_context": l2,
            "core_context": core_combined,
            "ensemble_baseline_core_context": ensemble,
        }
        kind_metrics = metrics_by_kind.setdefault(
            transformation,
            {method: _empty_metrics() for method in methods},
        )
        method_results = {}
        for method, candidates in candidate_sets.items():
            exact_rank, center_rank = _record(metrics[method], candidates, expected)
            _record(kind_metrics[method], candidates, expected)
            method_results[method] = {
                "candidate_count": len(candidates),
                "exact_precursor_rank": exact_rank,
                "center_rank": center_rank,
                "precursor_smiles": [
                    candidate.precursor_smiles for candidate in candidates
                ],
            }
        target_results.append(
            {
                "reaction_id": str(row.get("reaction_id") or ""),
                "reference_group": str(row.get("reference_id") or ""),
                "sampling_cohort": str(row.get("_sampling_cohort") or ""),
                "bond_kind": transformation,
                "target_smiles": target,
                "expected_precursor_smiles": expected,
                "expected_center_transition_key": transformation,
                "methods": method_results,
            }
        )
    report = {
        "definition_id": "diverse_core_retrosynthesis_benchmark.v1",
        "split": {
            "source_rows": len(rows),
            "training_rows": len(train),
            "held_out_rows": len(test),
            "evaluated_targets": len(selected_cases),
            "held_out_support_groups": len(
                {str(row.get("reference_id") or "") for row in test}
            ),
        },
        "builds": {
            "baseline": baseline_library.to_dict() | {"templates": None},
            "core": core_library.to_dict() | {"templates": None},
        },
        "eligibility_rejections": dict(sorted(eligibility.items())),
        "metrics": {method: _finalize(value) for method, value in metrics.items()},
        "metrics_by_bond": {
            kind: {
                method: _finalize(value)
                for method, value in method_values.items()
            }
            for kind, method_values in sorted(metrics_by_kind.items())
        },
        "target_results": target_results,
    }
    (output / "comparison.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
        newline="\n",
    )
    return report


__all__ = [
    "DIVERSE_COHORTS",
    "load_diverse_rows",
    "run_diverse_benchmark",
]
