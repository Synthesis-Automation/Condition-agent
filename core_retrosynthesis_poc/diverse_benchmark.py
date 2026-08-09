"""Leakage-safe benchmark across structurally diverse reaction archetypes."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from itertools import cycle
from pathlib import Path
from typing import Any, Dict, Sequence

from retrosynthesis_poc.library import iter_rows

from .comparison import split_by_reference
from .ensemble import merge_baseline_first
from .generic_compiler import (
    classify_reaction_with_site,
    compile_generic_templates,
)
from .generic_library import (
    build_generic_library,
    save_generic_library,
)
from .generic_search import disconnect_generic_target, rank_site_diverse


DIVERSE_COHORTS = {
    "acyl_substitution": ("Amide_formation*.jsonl.gz",),
    "c_c_coupling": ("suzuki_miyaura*.jsonl.gz",),
    "carbonyl_oxidation": ("DessMartin_periodinane*.jsonl.gz",),
    "carbonyl_reduction": ("NaBH4_carbonyl_reductions*.jsonl.gz",),
    "conjugate_addition": ("Michael_addition*.jsonl.gz",),
    "carbonyl_condensation": ("Reductive_amination*.jsonl.gz",),
    "ring_formation": ("CuAAC_azidealkyne*.jsonl.gz",),
}

STRESS_COHORTS = {
    "c_c_coupling": ("suzuki_miyaura*.jsonl.gz",),
    "conjugate_addition": ("Michael_addition*.jsonl.gz",),
}


def load_cohort_rows(
    source: str | Path,
    *,
    cohorts: Dict[str, Sequence[str]],
    max_rows_per_cohort: int = 200,
) -> tuple[Dict[str, Any], ...]:
    """Round-robin sample named cohorts without using names for routing."""

    if max_rows_per_cohort < 1:
        raise ValueError("max rows per cohort must be positive")
    values = []
    root = Path(source)
    for cohort, patterns in cohorts.items():
        files = tuple(
            sorted(
                {path for pattern in patterns for path in root.glob(pattern)},
                key=lambda path: path.as_posix(),
            )
        )
        iterators = {path: iter(iter_rows(path)) for path in files}
        selected_count = 0
        for path in cycle(files):
            if not iterators or selected_count >= max_rows_per_cohort:
                break
            iterator = iterators.get(path)
            if iterator is None:
                continue
            try:
                row = next(iterator)
            except StopIteration:
                del iterators[path]
                continue
            values.append({**row, "_sampling_cohort": cohort})
            selected_count += 1
    return tuple(values)


def load_diverse_rows(
    source: str | Path,
    *,
    max_rows_per_cohort: int = 200,
) -> tuple[Dict[str, Any], ...]:
    """Sample the seven diverse structural evaluation cohorts."""

    return load_cohort_rows(
        source,
        cohorts=DIVERSE_COHORTS,
        max_rows_per_cohort=max_rows_per_cohort,
    )


def load_stress_rows(
    source: str | Path,
    *,
    max_rows_per_cohort: int = 1_000,
) -> tuple[Dict[str, Any], ...]:
    """Sample larger C-C coupling and conjugate-addition stress cohorts."""

    return load_cohort_rows(
        source,
        cohorts=STRESS_COHORTS,
        max_rows_per_cohort=max_rows_per_cohort,
    )


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
        "top1_site": 0,
        "top5_site": 0,
        "top10_site": 0,
    }


def _record(
    metrics: Dict[str, Any],
    candidates: Sequence[Any],
    expected_precursors: str,
    expected_center: str,
    expected_site: str,
) -> tuple[int | None, int | None, int | None]:
    metrics["targets"] += 1
    metrics["targets_with_candidates"] += int(bool(candidates))
    metrics["candidate_count"] += len(candidates)
    metrics["valid_candidate_count"] += sum(
        candidate.forward_validation_status not in {"invalid", "unresolved"}
        for candidate in candidates
    )
    values = [candidate.precursor_smiles for candidate in candidates]
    exact_rank = values.index(expected_precursors) + 1 if expected_precursors in values else None
    center_values = [candidate.center_transition_key for candidate in candidates]
    center_rank = (
        center_values.index(expected_center) + 1
        if expected_center and expected_center in center_values
        else None
    )
    site_values = [candidate.disconnection_site_key for candidate in candidates]
    site_rank = (
        site_values.index(expected_site) + 1
        if expected_site and expected_site in site_values
        else None
    )
    for limit in (1, 5, 10):
        metrics[f"top{limit}_exact_precursor"] += int(
            exact_rank is not None and exact_rank <= limit
        )
        metrics[f"top{limit}_center"] += int(
            center_rank is not None and center_rank <= limit
        )
        metrics[f"top{limit}_site"] += int(
            site_rank is not None and site_rank <= limit
        )
    return exact_rank, center_rank, site_rank


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
        value[f"top{limit}_site_recall"] = metrics[f"top{limit}_site"] / targets
    return value


def _balanced_targets(
    cases: Sequence[tuple[Dict[str, Any], str, str, str, str, str]],
    max_targets_per_transformation: int,
) -> tuple[tuple[Dict[str, Any], str, str, str, str, str], ...]:
    grouped: dict[str, list[tuple[Dict[str, Any], str, str, str, str, str]]] = {}
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
    include_level_ablations: bool = True,
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
                str((row.get("reaction_core") or {}).get("center_transition_key") or ""),
                classify_reaction_with_site(precedent.mapped_reaction_smiles)[1],
            )
        )
    selected_cases = _balanced_targets(cases, max_targets_per_transformation)
    methods = ["baseline"]
    if include_level_ablations:
        methods.extend(("core_l1_context", "core_l2_context"))
    methods.extend(
        (
            "core_context",
            "core_site_diverse",
            "ensemble_baseline_core_context",
            "ensemble_baseline_core_site_diverse",
        )
    )
    methods = tuple(methods)
    metrics = {method: _empty_metrics() for method in methods}
    metrics_by_kind: dict[str, dict[str, Dict[str, Any]]] = {}
    metrics_by_difficulty: dict[str, dict[str, Dict[str, Any]]] = {}
    target_results = []
    failure_modes = {method: Counter() for method in methods}
    for (
        row,
        target,
        expected,
        transformation,
        expected_center,
        expected_site,
    ) in selected_cases:
        baseline = disconnect_generic_target(
            target,
            baseline_library,
            transformations=(transformation,),
            top_k=top_k,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        core_pool = disconnect_generic_target(
            target,
            core_library,
            transformations=(transformation,),
            levels=("L1", "L2"),
            top_k=max(top_k, max_candidates_to_validate),
            max_candidates_to_validate=max_candidates_to_validate,
        )
        core_combined = core_pool[:top_k]
        core_site_diverse = rank_site_diverse(core_pool)[:top_k]
        core_site_count = len(
            {
                candidate.disconnection_site_key
                for candidate in core_pool
                if candidate.disconnection_site_key
            }
        )
        difficulty = (
            "no_core_candidate"
            if not core_pool
            else "multi_site"
            if core_site_count > 1
            else "single_site"
        )
        ensemble = merge_baseline_first(baseline, core_combined, top_k)
        site_ensemble = merge_baseline_first(baseline, core_site_diverse, top_k)
        candidate_sets = {
            "baseline": baseline,
            "core_context": core_combined,
            "core_site_diverse": core_site_diverse,
            "ensemble_baseline_core_context": ensemble,
            "ensemble_baseline_core_site_diverse": site_ensemble,
        }
        if include_level_ablations:
            candidate_sets["core_l1_context"] = disconnect_generic_target(
                target,
                core_library,
                transformations=(transformation,),
                levels=("L1",),
                top_k=top_k,
                max_candidates_to_validate=max_candidates_to_validate,
            )
            candidate_sets["core_l2_context"] = disconnect_generic_target(
                target,
                core_library,
                transformations=(transformation,),
                levels=("L2",),
                top_k=top_k,
                max_candidates_to_validate=max_candidates_to_validate,
            )
        kind_metrics = metrics_by_kind.setdefault(
            transformation,
            {method: _empty_metrics() for method in methods},
        )
        difficulty_metrics = metrics_by_difficulty.setdefault(
            difficulty,
            {method: _empty_metrics() for method in methods},
        )
        method_results = {}
        for method, candidates in candidate_sets.items():
            exact_rank, center_rank, site_rank = _record(
                metrics[method],
                candidates,
                expected,
                expected_center,
                expected_site,
            )
            _record(
                kind_metrics[method],
                candidates,
                expected,
                expected_center,
                expected_site,
            )
            _record(
                difficulty_metrics[method],
                candidates,
                expected,
                expected_center,
                expected_site,
            )
            outcome = (
                "exact_top1"
                if exact_rank == 1
                else "exact_lower_rank"
                if exact_rank is not None
                else "correct_site_alternative_precursors"
                if site_rank is not None
                else "wrong_site_only"
                if candidates
                else "no_candidates"
            )
            failure_modes[method][outcome] += 1
            method_results[method] = {
                "candidate_count": len(candidates),
                "exact_precursor_rank": exact_rank,
                "center_rank": center_rank,
                "site_rank": site_rank,
                "outcome": outcome,
                "precursor_smiles": [
                    candidate.precursor_smiles for candidate in candidates
                ],
                "disconnection_site_keys": [
                    candidate.disconnection_site_key for candidate in candidates
                ],
            }
        target_results.append(
            {
                "reaction_id": str(row.get("reaction_id") or ""),
                "reference_group": str(row.get("reference_id") or ""),
                "sampling_cohort": str(row.get("_sampling_cohort") or ""),
                "bond_kind": transformation,
                "difficulty": difficulty,
                "core_candidate_site_count": core_site_count,
                "target_smiles": target,
                "expected_precursor_smiles": expected,
                "expected_center_transition_key": expected_center,
                "expected_disconnection_site_key": expected_site,
                "methods": method_results,
            }
        )
    report = {
        "definition_id": "diverse_core_retrosynthesis_benchmark.v2",
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
        "failure_modes": {
            method: dict(sorted(counts.items()))
            for method, counts in failure_modes.items()
        },
        "metrics_by_bond": {
            kind: {
                method: _finalize(value)
                for method, value in method_values.items()
            }
            for kind, method_values in sorted(metrics_by_kind.items())
        },
        "metrics_by_difficulty": {
            difficulty: {
                method: _finalize(value)
                for method, value in method_values.items()
            }
            for difficulty, method_values in sorted(metrics_by_difficulty.items())
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
    "STRESS_COHORTS",
    "load_cohort_rows",
    "load_diverse_rows",
    "load_stress_rows",
    "run_diverse_benchmark",
]
