"""Coverage-first benchmark for unrestricted data-derived graph operators."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import replace
from itertools import cycle
from pathlib import Path
from typing import Any, Dict, Sequence

from retrosynthesis_poc.library import iter_rows

from .comparison import split_by_reference
from .ensemble import merge_baseline_first
from .generic_compiler import analyze_generic_reaction, compile_generic_templates
from .generic_library import build_generic_library, save_generic_library
from .generic_models import GenericTemplateLibrary
from .generic_search import disconnect_generic_target


def _precedent_identity(precedent: Any) -> str:
    return str(precedent.reaction_id or precedent.mapped_reaction_smiles)


def load_operator_rows(
    source: str | Path,
    *,
    max_rows: int = 1_000,
) -> tuple[Dict[str, Any], ...]:
    """Round-robin sample every dataset shard without family-based routing."""

    if max_rows < 1:
        raise ValueError("max rows must be positive")
    root = Path(source)
    files = (
        (root,)
        if root.is_file()
        else tuple(
            sorted(
                root.glob("*.jsonl.gz"),
                key=lambda path: hashlib.sha256(path.name.encode()).digest(),
            )
        )
    )
    iterators = {path: iter(iter_rows(path)) for path in files}
    values = []
    for path in cycle(files):
        if not iterators or len(values) >= max_rows:
            break
        iterator = iterators.get(path)
        if iterator is None:
            continue
        try:
            row = next(iterator)
        except StopIteration:
            del iterators[path]
            continue
        values.append({**row, "_sampling_shard": path.name})
    return tuple(values)


def _supported_view(library: GenericTemplateLibrary) -> GenericTemplateLibrary:
    templates = tuple(
        template
        for template in library.templates
        if template.transformation_kind is not None
        and template.abstraction_level in {"L1", "L2"}
    )
    operator_ids = {template.operator_id for template in templates}
    accepted_ids = {
        _precedent_identity(precedent)
        for template in templates
        for precedent in template.precedents
    }
    rejections = Counter(library.rejection_counts)
    rejections["unsupported_edit_archetype"] += max(
        0,
        library.accepted_observation_count - len(accepted_ids),
    )
    return replace(
        library,
        templates=templates,
        operators=tuple(
            operator
            for operator in library.operators
            if operator.operator_id in operator_ids
        ),
        accepted_observation_count=len(accepted_ids),
        rejection_counts=dict(sorted(rejections.items())),
        definition={
            **library.definition,
            "admission_mode": "supported_view",
            "levels": ["L1", "L2"],
        },
    )


def _census(library: GenericTemplateLibrary) -> Dict[str, Any]:
    annotated_reactions = {
        _precedent_identity(precedent)
        for template in library.templates
        if template.named_annotations
        for precedent in template.precedents
    }
    unannotated_reactions = {
        _precedent_identity(precedent)
        for template in library.templates
        if not template.named_annotations
        for precedent in template.precedents
    }
    annotated_count = int(
        library.definition["annotated_accepted_observation_count"]
        if "annotated_accepted_observation_count" in library.definition
        else len(annotated_reactions)
    )
    unannotated_count = int(
        library.definition["unannotated_accepted_observation_count"]
        if "unannotated_accepted_observation_count" in library.definition
        else len(unannotated_reactions)
    )
    operator_rows = sorted(
        (
            {
                "operator_id": operator.operator_id,
                "operator_signature": operator.operator_signature,
                "observation_support": operator.observation_support,
                "independent_reference_support": (
                    operator.independent_reference_support
                ),
                "realization_count": len(operator.realization_ids),
                "levels": list(operator.abstraction_levels),
                "named_annotations": list(operator.named_annotations),
            }
            for operator in library.operators
        ),
        key=lambda value: (
            -int(value["independent_reference_support"]),
            str(value["operator_id"]),
        ),
    )
    return {
        "source_rows": library.source_row_count,
        "accepted_observations": library.accepted_observation_count,
        "acceptance_fraction": (
            library.accepted_observation_count / library.source_row_count
            if library.source_row_count
            else 0.0
        ),
        "template_count": len(library.templates),
        "operator_count": len(library.operators),
        "realization_count": len(
            {
                template.realization_id
                for template in library.templates
                if template.realization_id
            }
        ),
        "annotated_observation_count": annotated_count,
        "unannotated_observation_count": unannotated_count,
        "retained_precedent_observation_count": len(
            annotated_reactions.union(unannotated_reactions)
        ),
        "rejection_counts": dict(sorted(library.rejection_counts.items())),
        "operators": operator_rows,
    }


def _empty_metrics() -> Dict[str, Any]:
    value: Dict[str, Any] = {
        "targets": 0,
        "targets_with_candidates": 0,
        "candidate_count": 0,
        "valid_candidate_count": 0,
        "supported_alternative_targets": 0,
    }
    for label in ("exact_precursor", "synthon", "operator", "site"):
        for limit in (1, 5, 10, 25):
            value[f"top{limit}_{label}"] = 0
    return value


def _rank(values: Sequence[str], expected: str) -> int | None:
    return values.index(expected) + 1 if expected and expected in values else None


def _record(
    metrics: Dict[str, Any],
    candidates: Sequence[Any],
    *,
    expected_precursors: str,
    expected_synthon: str,
    expected_operator: str,
    expected_site: str,
) -> Dict[str, Any]:
    metrics["targets"] += 1
    metrics["targets_with_candidates"] += int(bool(candidates))
    metrics["candidate_count"] += len(candidates)
    metrics["valid_candidate_count"] += sum(
        candidate.forward_validation_status not in {"invalid", "unresolved"}
        for candidate in candidates
    )
    exact_rank = _rank(
        [candidate.precursor_smiles for candidate in candidates],
        expected_precursors,
    )
    synthon_rank = _rank(
        [candidate.synthon_signature for candidate in candidates],
        expected_synthon,
    )
    operator_rank = _rank(
        [candidate.operator_id for candidate in candidates],
        expected_operator,
    )
    site_rank = _rank(
        [candidate.disconnection_site_key for candidate in candidates],
        expected_site,
    )
    ranks = {
        "exact_precursor": exact_rank,
        "synthon": synthon_rank,
        "operator": operator_rank,
        "site": site_rank,
    }
    for label, rank in ranks.items():
        for limit in (1, 5, 10, 25):
            metrics[f"top{limit}_{label}"] += int(
                rank is not None and rank <= limit
            )
    supported_alternative = any(
        candidate.precursor_smiles != expected_precursors
        and candidate.disconnection_site_key == expected_site
        and candidate.operator_id == expected_operator
        and candidate.synthon_signature == expected_synthon
        for candidate in candidates
    )
    metrics["supported_alternative_targets"] += int(supported_alternative)
    correct_site_operator = any(
        candidate.disconnection_site_key == expected_site
        and candidate.operator_id == expected_operator
        for candidate in candidates
    )
    outcome = (
        "exact_top1"
        if exact_rank == 1
        else "exact_lower_rank"
        if exact_rank is not None
        else "synthon_equivalent_alternative"
        if supported_alternative
        else "correct_site_and_operator"
        if correct_site_operator
        else "correct_site_other_operator"
        if site_rank is not None
        else "correct_operator_other_site"
        if operator_rank is not None
        else "other_valid_candidates"
        if candidates
        else "no_candidates"
    )
    return {
        "exact_precursor_rank": exact_rank,
        "synthon_rank": synthon_rank,
        "operator_rank": operator_rank,
        "site_rank": site_rank,
        "supported_alternative": supported_alternative,
        "outcome": outcome,
    }


def _finalize(metrics: Dict[str, Any]) -> Dict[str, Any]:
    targets = max(1, int(metrics["targets"]))
    candidates = max(1, int(metrics["candidate_count"]))
    value = dict(metrics)
    value["target_coverage"] = metrics["targets_with_candidates"] / targets
    value["mean_candidates_per_target"] = metrics["candidate_count"] / targets
    value["valid_candidate_fraction"] = (
        metrics["valid_candidate_count"] / candidates
    )
    value["supported_alternative_rate"] = (
        metrics["supported_alternative_targets"] / targets
    )
    for label in ("exact_precursor", "synthon", "operator", "site"):
        for limit in (1, 5, 10, 25):
            value[f"top{limit}_{label}_recall"] = (
                metrics[f"top{limit}_{label}"] / targets
            )
    return value


def _select_cases(
    cases: Sequence[Dict[str, Any]],
    *,
    max_targets: int,
    max_targets_per_operator: int,
) -> tuple[Dict[str, Any], ...]:
    by_operator: dict[str, list[Dict[str, Any]]] = {}
    for case in cases:
        by_operator.setdefault(case["expected_operator_id"], []).append(case)
    selected = []
    for operator, values in sorted(by_operator.items()):
        ranked = sorted(
            values,
            key=lambda value: hashlib.sha256(
                f"{operator}\0{value['reaction_id']}".encode()
            ).digest(),
        )
        selected.extend(ranked[:max_targets_per_operator])
    annotated = [case for case in selected if case["named_annotation"]]
    unannotated = [case for case in selected if not case["named_annotation"]]
    for values in (annotated, unannotated):
        values.sort(
            key=lambda value: hashlib.sha256(
                f"operator-target\0{value['reaction_id']}".encode()
            ).digest()
        )
    balanced = []
    for index in range(max(len(annotated), len(unannotated))):
        for values in (unannotated, annotated):
            if index < len(values):
                balanced.append(values[index])
                if len(balanced) >= max_targets:
                    return tuple(balanced)
    return tuple(balanced)


def _method_payload(candidates: Sequence[Any], ranks: Dict[str, Any]) -> Dict[str, Any]:
    return {
        **ranks,
        "candidate_count": len(candidates),
        "precursor_smiles": [candidate.precursor_smiles for candidate in candidates],
        "disconnection_site_keys": [
            candidate.disconnection_site_key for candidate in candidates
        ],
        "operator_signatures": [
            candidate.operator_signature for candidate in candidates
        ],
        "operator_ids": [candidate.operator_id for candidate in candidates],
        "synthon_signatures": [
            candidate.synthon_signature for candidate in candidates
        ],
        "retrieval_levels": [
            candidate.abstraction_level for candidate in candidates
        ],
    }


def run_operator_coverage_benchmark(
    rows: Sequence[Dict[str, Any]],
    output_directory: str | Path,
    *,
    test_fraction: float = 0.25,
    max_targets: int = 80,
    max_targets_per_operator: int = 3,
    top_k: int = 25,
    max_templates_to_apply: int = 500,
    max_candidates_to_validate: int = 100,
) -> Dict[str, Any]:
    """Build unrestricted operators and evaluate coverage-equivalence levels."""

    train, test = split_by_reference(
        rows,
        test_fraction=test_fraction,
        seed="data-derived-operator-poc-v1",
    )
    library = build_generic_library(
        train,
        engine="reaction_core",
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
    )
    supported = _supported_view(library)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    save_generic_library(library, output / "operator_library.json.gz")
    save_generic_library(supported, output / "supported_library.json.gz")
    census = _census(library)
    (output / "operator_census.json").write_text(
        json.dumps(census, indent=2, sort_keys=True),
        encoding="utf-8",
        newline="\n",
    )

    cases = []
    eligibility = Counter()
    for row in test:
        result = compile_generic_templates(
            row,
            engine="reaction_core",
            levels=("L0", "L1", "L2"),
            admission_mode="data_driven",
        )
        if not result.templates:
            eligibility[str(result.rejection_reason or "unknown_rejection")] += 1
            continue
        template = result.templates[0]
        precedent = template.precedents[0]
        identity = analyze_generic_reaction(precedent.mapped_reaction_smiles)
        if identity is None or not identity.disconnection_site_key:
            eligibility["expected_identity_unavailable"] += 1
            continue
        cases.append(
            {
                "row": row,
                "reaction_id": precedent.reaction_id,
                "reference_group": precedent.reference_id,
                "target_smiles": precedent.product_smiles,
                "expected_precursor_smiles": precedent.precursor_smiles,
                "expected_site": identity.disconnection_site_key,
                "expected_operator_signature": template.operator_signature,
                "expected_operator_id": template.operator_id,
                "expected_synthon": template.synthon_signature,
                "named_annotation": template.transformation_kind,
            }
        )
    selected = _select_cases(
        cases,
        max_targets=max_targets,
        max_targets_per_operator=max_targets_per_operator,
    )
    methods = (
        "supported_l1_l2",
        "data_l2",
        "data_l2_l1",
        "data_ladder",
    )
    metrics = {method: _empty_metrics() for method in methods}
    metrics_by_scope = {
        scope: {method: _empty_metrics() for method in methods}
        for scope in ("annotated", "unannotated")
    }
    outcomes = {method: Counter() for method in methods}
    target_results = []
    for case in selected:
        target = case["target_smiles"]
        data_l2 = disconnect_generic_target(
            target,
            library,
            levels=("L2",),
            top_k=top_k,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        data_l1 = disconnect_generic_target(
            target,
            library,
            levels=("L1",),
            top_k=top_k,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        data_l0 = disconnect_generic_target(
            target,
            library,
            levels=("L0",),
            top_k=top_k,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
        )
        data_l2_l1 = merge_baseline_first(data_l2, data_l1, top_k)
        data_ladder = merge_baseline_first(data_l2_l1, data_l0, top_k)
        candidate_sets = {
            "supported_l1_l2": disconnect_generic_target(
                target,
                supported,
                levels=("L1", "L2"),
                top_k=top_k,
                max_templates_to_apply=max_templates_to_apply,
                max_candidates_to_validate=max_candidates_to_validate,
            ),
            "data_l2": data_l2,
            "data_l2_l1": data_l2_l1,
            "data_ladder": data_ladder,
        }
        scope = "annotated" if case["named_annotation"] else "unannotated"
        method_results = {}
        for method, candidates in candidate_sets.items():
            record_args = {
                "expected_precursors": case["expected_precursor_smiles"],
                "expected_synthon": case["expected_synthon"],
                "expected_operator": case["expected_operator_id"],
                "expected_site": case["expected_site"],
            }
            ranks = _record(metrics[method], candidates, **record_args)
            _record(metrics_by_scope[scope][method], candidates, **record_args)
            outcomes[method][str(ranks["outcome"])] += 1
            method_results[method] = _method_payload(candidates, ranks)
        target_results.append(
            {
                "reaction_id": case["reaction_id"],
                "reference_group": case["reference_group"],
                "bond_kind": case["named_annotation"] or "unannotated_operator",
                "operator_scope": scope,
                "target_smiles": target,
                "expected_precursor_smiles": case["expected_precursor_smiles"],
                "expected_disconnection_site_key": case["expected_site"],
                "expected_operator_signature": case[
                    "expected_operator_signature"
                ],
                "expected_operator_id": case["expected_operator_id"],
                "expected_synthon_signature": case["expected_synthon"],
                "condition_compatibility_status": "not_evaluated_no_query_conditions",
                "methods": method_results,
            }
        )
    report = {
        "definition_id": "data_derived_operator_coverage_benchmark.v1",
        "split": {
            "source_rows": len(rows),
            "training_rows": len(train),
            "held_out_rows": len(test),
            "held_out_support_groups": len(
                {
                    str(row.get("reference_id") or row.get("reaction_id") or "")
                    for row in test
                }
            ),
            "eligible_held_out_targets": len(cases),
            "evaluated_targets": len(selected),
            "evaluated_annotated_targets": sum(
                bool(case["named_annotation"]) for case in selected
            ),
            "evaluated_unannotated_targets": sum(
                not case["named_annotation"] for case in selected
            ),
        },
        "census": census,
        "supported_view": _census(supported),
        "eligibility_rejections": dict(sorted(eligibility.items())),
        "equivalence_levels": {
            "site": "same edited product atoms and bonds",
            "operator": "same handle-independent normalized graph edits",
            "synthon": "same mapped precursor skeletons after handle removal",
            "exact": "same canonical recorded precursor molecules",
            "supported_alternative": (
                "same site, operator, and synthon with a different source-supported "
                "realization"
            ),
            "condition_compatibility": "not evaluated without query conditions",
        },
        "metrics": {method: _finalize(value) for method, value in metrics.items()},
        "metrics_by_scope": {
            scope: {
                method: _finalize(value)
                for method, value in method_values.items()
            }
            for scope, method_values in metrics_by_scope.items()
        },
        "outcomes": {
            method: dict(sorted(counts.items()))
            for method, counts in outcomes.items()
        },
        "target_results": target_results,
    }
    (output / "coverage_comparison.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
        newline="\n",
    )
    return report


__all__ = ["load_operator_rows", "run_operator_coverage_benchmark"]
