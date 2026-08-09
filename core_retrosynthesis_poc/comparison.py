"""Leakage-safe paired evaluation of baseline and core-derived templates."""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Sequence

from reactive_taxonomy import featurize_reaction
from retrosynthesis_poc import build_library as build_baseline_library
from retrosynthesis_poc import disconnect_target as baseline_disconnect
from retrosynthesis_poc import save_library as save_baseline_library
from retrosynthesis_poc.extraction import extract_cx_template

from .compiler import compile_core_templates
from .ensemble import merge_baseline_first
from .library import build_library as build_core_library
from .library import save_library as save_core_library
from .search import disconnect_target as core_disconnect


def _support_group(row: Dict[str, Any]) -> str:
    reference = str(row.get("reference_id") or "").strip()
    if reference:
        return f"reference:{reference}"
    identity = row.get("reference_identity") or {}
    nested = str(identity.get("reference_id") or "").strip()
    if nested:
        return f"reference:{nested}"
    core = row.get("reaction_core") or {}
    return str(
        core.get("mapping_equivalence_key")
        or core.get("core_id")
        or row.get("reaction_id")
        or row.get("observation_id")
        or "unidentified"
    )


def split_by_reference(
    rows: Sequence[Dict[str, Any]],
    *,
    test_fraction: float = 0.2,
    seed: str = "core-retro-poc-v1",
) -> tuple[tuple[Dict[str, Any], ...], tuple[Dict[str, Any], ...]]:
    """Make a deterministic reference-group split with no group leakage."""

    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test fraction must be between zero and one")
    grouped: dict[str, list[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(_support_group(row), []).append(row)
    ranked_groups = sorted(
        grouped,
        key=lambda group: hashlib.sha256(
            f"{seed}\0{group}".encode("utf-8")
        ).digest(),
    )
    if len(ranked_groups) < 2:
        return tuple(rows), ()
    test_group_count = max(
        1,
        min(len(ranked_groups) - 1, round(len(ranked_groups) * test_fraction)),
    )
    test_groups = set(ranked_groups[:test_group_count])
    train = []
    test = []
    for group, group_rows in grouped.items():
        (test if group in test_groups else train).extend(group_rows)
    return tuple(train), tuple(test)


def _expected_case(row: Dict[str, Any]) -> tuple[str, str, str] | None:
    baseline = extract_cx_template(row)
    if baseline.template is not None:
        precedent = baseline.template.precedent
        return (
            precedent.product_smiles,
            precedent.precursor_smiles,
            baseline.template.bond_kind,
        )
    core = compile_core_templates(row, levels=("L1",))
    if core.templates:
        template = core.templates[0]
        precedent = template.precedents[0]
        return precedent.product_smiles, precedent.precursor_smiles, template.bond_kind
    return None


@lru_cache(maxsize=50_000)
def _center_key(reaction_smiles: str) -> str:
    analysis = featurize_reaction(reaction_smiles)
    return analysis.reaction_core.center_transition_key if analysis.reaction_core else ""


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
    expected_center: str,
) -> None:
    metrics["targets"] += 1
    metrics["targets_with_candidates"] += int(bool(candidates))
    metrics["candidate_count"] += len(candidates)
    metrics["valid_candidate_count"] += sum(
        candidate.forward_validation_status not in {"invalid", "unresolved"}
        for candidate in candidates
    )
    precursor_values = [candidate.precursor_smiles for candidate in candidates]
    center_values = [
        str(getattr(candidate, "center_transition_key", ""))
        or _center_key(candidate.proposed_reaction_smiles)
        for candidate in candidates
    ]
    for limit in (1, 5, 10):
        metrics[f"top{limit}_exact_precursor"] += int(
            expected_precursors in precursor_values[:limit]
        )
        metrics[f"top{limit}_center"] += int(
            bool(expected_center) and expected_center in center_values[:limit]
        )


def _finalize(metrics: Dict[str, Any]) -> Dict[str, Any]:
    targets = max(1, int(metrics["targets"]))
    candidates = max(1, int(metrics["candidate_count"]))
    result = dict(metrics)
    result["target_coverage"] = metrics["targets_with_candidates"] / targets
    result["mean_candidates_per_target"] = metrics["candidate_count"] / targets
    result["valid_candidate_fraction"] = metrics["valid_candidate_count"] / candidates
    for limit in (1, 5, 10):
        result[f"top{limit}_exact_precursor_recall"] = (
            metrics[f"top{limit}_exact_precursor"] / targets
        )
        result[f"top{limit}_center_recall"] = metrics[f"top{limit}_center"] / targets
    return result


def _balanced_evaluation_rows(
    rows: Sequence[Dict[str, Any]],
    cases_by_identity: Dict[int, tuple[str, str, str]],
    limit: int | None,
) -> tuple[Dict[str, Any], ...]:
    if limit is None:
        return tuple(rows)
    buckets = {
        bond: [row for row in rows if cases_by_identity[id(row)][2] == bond]
        for bond in ("C-N", "C-O", "C-S")
    }
    selected = []
    while len(selected) < limit and any(buckets.values()):
        for bond in ("C-N", "C-O", "C-S"):
            if buckets[bond] and len(selected) < limit:
                selected.append(buckets[bond].pop(0))
    return tuple(selected)


def _rank(values: Sequence[str], expected: str) -> int | None:
    try:
        return values.index(expected) + 1
    except ValueError:
        return None


def run_comparison(
    rows: Iterable[Dict[str, Any]],
    output_directory: str | Path,
    *,
    test_fraction: float = 0.2,
    top_k: int = 10,
    max_test_targets: int | None = None,
    max_candidates_to_validate: int = 20,
) -> Dict[str, Any]:
    """Build both methods on training references and evaluate held-out targets."""

    if max_test_targets is not None and max_test_targets < 1:
        raise ValueError("max test targets must be positive")
    if max_candidates_to_validate < top_k:
        raise ValueError("validation shortlist must be at least top-k")
    values = tuple(rows)
    eligible_cases = tuple(
        (row, case)
        for row in values
        if (case := _expected_case(row)) is not None
    )
    eligible_rows = tuple(row for row, _case in eligible_cases)
    eligible_train, test = split_by_reference(
        eligible_rows,
        test_fraction=test_fraction,
    )
    if not eligible_train or not test:
        raise ValueError(
            "comparison requires at least two eligible reference groups"
        )
    test_groups = {_support_group(row) for row in test}
    train = tuple(
        row for row in values if _support_group(row) not in test_groups
    )
    train_groups = {_support_group(row) for row in train}
    if train_groups.intersection(test_groups):
        raise RuntimeError("reference split leaked between train and test")
    baseline_library, baseline_build = build_baseline_library(train)
    core_library, core_build = build_core_library(train)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    save_baseline_library(baseline_library, output / "baseline_templates.json.gz")
    save_core_library(core_library, output / "core_templates.json.gz")

    metrics = {
        "baseline": _empty_metrics(),
        "ensemble_baseline_l1_context": _empty_metrics(),
        "core_l1_neutral": _empty_metrics(),
        "core_l1_context": _empty_metrics(),
        "core_l2_neutral": _empty_metrics(),
        "core_l2_context": _empty_metrics(),
        "core_neutral": _empty_metrics(),
        "core_context": _empty_metrics(),
    }
    metrics_by_bond = {
        method: {bond: _empty_metrics() for bond in ("C-N", "C-O", "C-S")}
        for method in metrics
    }
    cases_by_identity = {id(row): case for row, case in eligible_cases}
    evaluated = 0
    target_results = []
    evaluation_rows = _balanced_evaluation_rows(
        test,
        cases_by_identity,
        max_test_targets,
    )
    for row in evaluation_rows:
        case = cases_by_identity[id(row)]
        target, expected_precursors, bond_kind = case
        expected_center = str(
            (row.get("reaction_core") or {}).get("center_transition_key") or ""
        )
        baseline_candidates = baseline_disconnect(
            target,
            baseline_library,
            allowed_bonds=(bond_kind,),
            top_k=top_k,
        )
        method_candidates = {"baseline": baseline_candidates}
        for level_name, levels in (
            ("l1", ("L1",)),
            ("l2", ("L2",)),
            ("", ("L1", "L2")),
        ):
            method_stem = f"core_{level_name}_" if level_name else "core_"
            for ranking_name, use_context in (
                ("neutral", False),
                ("context", True),
            ):
                method_candidates[f"{method_stem}{ranking_name}"] = (
                    core_disconnect(
                        target,
                        core_library,
                        allowed_bonds=(bond_kind,),
                        levels=levels,
                        top_k=top_k,
                        max_candidates_to_validate=max_candidates_to_validate,
                        use_context=use_context,
                    )
                )
        method_candidates["ensemble_baseline_l1_context"] = (
            merge_baseline_first(
                method_candidates["baseline"],
                method_candidates["core_l1_context"],
                top_k,
            )
        )
        for method, candidates in method_candidates.items():
            _record(
                metrics[method],
                candidates,
                expected_precursors,
                expected_center,
            )
            _record(
                metrics_by_bond[method][bond_kind],
                candidates,
                expected_precursors,
                expected_center,
            )
        target_results.append(
            {
                "reaction_id": str(row.get("reaction_id") or ""),
                "reference_group": _support_group(row),
                "bond_kind": bond_kind,
                "target_smiles": target,
                "expected_precursor_smiles": expected_precursors,
                "expected_center_transition_key": expected_center,
                "methods": {
                    method: {
                        "candidate_count": len(candidates),
                        "exact_precursor_rank": _rank(
                            [
                                candidate.precursor_smiles
                                for candidate in candidates
                            ],
                            expected_precursors,
                        ),
                        "center_rank": _rank(
                            [
                                str(
                                    getattr(
                                        candidate,
                                        "center_transition_key",
                                        "",
                                    )
                                )
                                or _center_key(
                                    candidate.proposed_reaction_smiles
                                )
                                for candidate in candidates
                            ],
                            expected_center,
                        ),
                        "precursor_smiles": [
                            candidate.precursor_smiles
                            for candidate in candidates
                        ],
                    }
                    for method, candidates in method_candidates.items()
                },
            }
        )
        evaluated += 1

    report = {
        "definition_id": "core_retrosynthesis_comparison.v1.2",
        "split": {
            "source_rows": len(values),
            "eligible_rows": len(eligible_rows),
            "training_rows": len(train),
            "held_out_rows": len(test),
            "training_support_groups": len(train_groups),
            "held_out_support_groups": len(test_groups),
            "evaluated_targets": evaluated,
            "max_test_targets": max_test_targets,
            "max_candidates_to_validate": max_candidates_to_validate,
            "test_fraction": test_fraction,
            "reference_leakage_count": 0,
        },
        "builds": {
            "baseline": {
                "accepted_observations": baseline_build.accepted_observation_count,
                "unique_templates": baseline_build.unique_template_count,
                "templates_by_bond": {
                    bond: sum(
                        template.bond_kind == bond
                        for template in baseline_library.templates
                    )
                    for bond in ("C-N", "C-O", "C-S")
                },
                "rejections": baseline_build.rejection_counts,
            },
            "core": {
                "accepted_observations": core_build.accepted_observation_count,
                "unique_templates": core_build.unique_template_count,
                "unique_operators": core_build.unique_operator_count,
                "templates_by_level": {
                    level: sum(
                        template.abstraction_level == level
                        for template in core_library.templates
                    )
                    for level in ("L1", "L2")
                },
                "operators_by_level": {
                    level: len(
                        {
                            template.operator_id
                            for template in core_library.templates
                            if template.abstraction_level == level
                        }
                    )
                    for level in ("L1", "L2")
                },
                "rejections": core_build.rejection_counts,
            },
        },
        "metrics": {name: _finalize(value) for name, value in metrics.items()},
        "metrics_by_bond": {
            method: {
                bond: _finalize(values)
                for bond, values in bond_metrics.items()
            }
            for method, bond_metrics in metrics_by_bond.items()
        },
        "target_results": target_results,
    }
    (output / "comparison.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = ["run_comparison", "split_by_reference"]
