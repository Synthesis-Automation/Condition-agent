"""Leakage-safe calibration and chemist review for reaction-core retrieval."""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from reactive_taxonomy.reaction_core.quality import (
    load_reaction_core_quality_rules,
)

from .core_retrieval import (
    load_reaction_core_retrieval_rules,
    retrieve_core_pool_with_trace,
)
from .evaluation import grouped_holdout_split
from .fallback_similarity import (
    assess_fallback_similarity,
    compatibility_signature_from_fallback,
    load_fallback_retrieval_rules,
)
from .generic_indexing import (
    GenericIndexedReaction,
    build_generic_index_from_rows,
    load_generic_index,
)
from .recipe_ranking import (
    load_generic_ranking_rules,
    rank_condition_recipes,
)
from .support import (
    load_evidence_support_rules,
    mapping_equivalence_key,
    summarize_evidence_support,
)


REACTION_CORE_EVALUATION_SCHEMA_VERSION = "1.0"
REACTION_CORE_EVALUATOR_VERSION = "reaction_core_calibration.v1"
_CORE_LEVELS = (
    "reaction_core_exact",
    "reaction_core_local",
    "reaction_core_context",
)
_REVIEW_FIELDS = (
    "case_id",
    "query_reaction_smiles",
    "query_reaction_label",
    "query_quality_status",
    "retrieval_level",
    "independent_compatible_candidate_count",
    "precedent_reaction_smiles",
    "precedent_reaction_labels",
    "chemist_core_match",
    "chemist_notes",
)


def _quality_status(core: Mapping[str, Any]) -> str:
    quality = core.get("quality")
    if isinstance(quality, Mapping):
        return str(quality.get("status") or "missing")
    return "missing"


def _reaction_label(value: Mapping[str, Any]) -> str:
    """Return the one terminal reaction label from a serialized contract."""
    return str(value.get("text") or "")


def _mean(values: Iterable[int]) -> float | None:
    numbers = tuple(values)
    return round(sum(numbers) / len(numbers), 6) if numbers else None


def _mapping_summary(rows: Iterable[GenericIndexedReaction]) -> dict[str, Any]:
    values = tuple(rows)
    groups: dict[str, list[GenericIndexedReaction]] = {}
    for row in values:
        key = mapping_equivalence_key(row)
        if key:
            groups.setdefault(key, []).append(row)
    support = summarize_evidence_support(values)
    return {
        "row_count": len(values),
        "rows_with_mapping_equivalence": sum(len(group) for group in groups.values()),
        "rows_without_mapping_equivalence": sum(
            not mapping_equivalence_key(row) for row in values
        ),
        "mapping_equivalence_group_count": len(groups),
        "multirow_mapping_equivalence_group_count": sum(
            len(group) > 1 for group in groups.values()
        ),
        "largest_mapping_equivalence_group_size": max(
            (len(group) for group in groups.values()),
            default=0,
        ),
        "cross_canonical_reaction_group_count": sum(
            len(
                {
                    row.canonical_reaction_id
                    for row in group
                    if row.canonical_reaction_id
                }
            )
            > 1
            for group in groups.values()
        ),
        "independent_support_count": support.independent_count,
        "mapping_deduplicated_independent_count": (
            support.mapping_deduplicated_independent_count
        ),
    }


def _stratified_metrics(
    cases: Iterable[Mapping[str, Any]],
    field: str,
) -> dict[str, dict[str, Any]]:
    groups: dict[str, list[Mapping[str, Any]]] = {}
    for case in cases:
        groups.setdefault(str(case.get(field) or "unknown"), []).append(case)
    result = {}
    for key, members in sorted(groups.items()):
        count = len(members)
        covered = sum(bool(member["valid"]) for member in members)
        candidate_count = sum(int(member["candidate_count"]) for member in members)
        compatible_count = sum(
            int(member["compatible_candidate_count"]) for member in members
        )
        result[key] = {
            "query_count": count,
            "covered_query_count": covered,
            "coverage_rate": round(covered / count, 6) if count else 0.0,
            "top1_recipe_recovery_rate": round(
                sum(bool(member["top1_recipe_match"]) for member in members)
                / count,
                6,
            )
            if count
            else 0.0,
            "topk_recipe_recovery_rate": round(
                sum(bool(member["topk_recipe_match"]) for member in members)
                / count,
                6,
            )
            if count
            else 0.0,
            "candidate_count": candidate_count,
            "compatible_candidate_count": compatible_count,
            "compatibility_survival_rate": round(
                compatible_count / candidate_count,
                6,
            )
            if candidate_count
            else 0.0,
            "mean_independent_compatible_support": _mean(
                int(member["independent_compatible_candidate_count"])
                for member in members
                if member["valid"]
            ),
        }
    return result


def _markdown(report: Mapping[str, Any]) -> str:
    metrics = report["metrics"]
    mapping = report["mapping_equivalence"]["all_rows"]
    lines = [
        "# Reaction-core calibration report",
        "",
        f"- Eligible held-out cores: {metrics['core_query_count']}",
        f"- Core retrieval coverage: {metrics['coverage_rate']:.2%}",
        f"- Top-1 recipe recovery: {metrics['top1_recipe_recovery_rate']:.2%}",
        f"- Top-k recipe recovery: {metrics['topk_recipe_recovery_rate']:.2%}",
        f"- Compatibility survival: {metrics['compatibility_survival_rate']:.2%}",
        f"- Mapping-equivalence groups: {mapping['mapping_equivalence_group_count']}",
        f"- Deduplicated support units: {mapping['mapping_deduplicated_independent_count']}",
        "",
        "## Retrieval tiers",
        "",
        "| Tier | Queries | Coverage | Top-1 | Top-k |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    by_level = report["by_retrieval_level"]
    for level in _CORE_LEVELS:
        value = by_level.get(level, {})
        lines.append(
            f"| {level} | {value.get('query_count', 0)} | "
            f"{value.get('coverage_rate', 0.0):.2%} | "
            f"{value.get('top1_recipe_recovery_rate', 0.0):.2%} | "
            f"{value.get('topk_recipe_recovery_rate', 0.0):.2%} |"
        )
    lines.extend(
        (
            "",
            "The accompanying chemist-review CSV omits expected recipes and",
            "contains blank adjudication columns to reduce review bias.",
            "",
        )
    )
    return "\n".join(lines)


def evaluate_reaction_core_index(
    records_path: str | Path,
    output_dir: str | Path,
    *,
    test_fraction: float = 0.2,
    seed: int = 17,
    top_k: int = 5,
    minimum_pool_size: int | None = None,
    split_mode: str = "grouped_random",
) -> dict[str, Any]:
    """Evaluate core-only retrieval with leakage-safe held-out precedents."""
    if top_k < 1:
        raise ValueError("top_k must be positive")
    full_index = load_generic_index(records_path)
    split = grouped_holdout_split(
        full_index.rows,
        test_fraction=test_fraction,
        seed=seed,
        split_mode=split_mode,  # type: ignore[arg-type]
    )
    train_index = build_generic_index_from_rows(split.train_rows)
    train_recipe_core_ids = {row.recipe_core_id for row in split.train_rows}
    cases = []
    review_rows = []
    unavailable = 0
    for row in split.test_rows:
        core = row.reaction_core
        if not core:
            unavailable += 1
            continue
        fallback = row.fallback_descriptor
        retrieval = retrieve_core_pool_with_trace(
            core,
            compatibility_signature_from_fallback(fallback),
            train_index,
            minimum_pool_size=minimum_pool_size,
        )
        recommendations = ()
        if retrieval.pool:
            query = {"reaction_core": core, "fallback_descriptor": fallback}

            def similarity(
                query_value: Mapping[str, Any],
                precedent: GenericIndexedReaction,
            ) -> Any:
                return assess_fallback_similarity(
                    query_value["fallback_descriptor"],
                    precedent.fallback_descriptor,
                )

            recommendations = rank_condition_recipes(
                query,
                retrieval.pool,
                retrieval_level=retrieval.level,
                top_k=top_k,
                similarity_assessor=similarity,
            )
        recommended_core_ids = tuple(
            recommendation.recipe_core_id for recommendation in recommendations
        )
        top1_match = bool(
            recommended_core_ids and recommended_core_ids[0] == row.recipe_core_id
        )
        topk_match = row.recipe_core_id in recommended_core_ids
        base_level = next(
            (
                level
                for level in _CORE_LEVELS
                if retrieval.level.startswith(level)
            ),
            retrieval.level,
        )
        case = {
            "canonical_reaction_id": row.canonical_reaction_id,
            "reaction_id": row.reaction_id,
            "observation_id": row.observation_id,
            "mapping_equivalence_key": mapping_equivalence_key(row),
            "quality_status": _quality_status(core),
            "actual_recipe_core_id": row.recipe_core_id,
            "recipe_seen_in_training": row.recipe_core_id in train_recipe_core_ids,
            "valid": bool(recommendations),
            "retrieval_level": base_level,
            "retrieval_result_level": retrieval.level,
            "candidate_count": retrieval.candidate_count,
            "independent_candidate_count": retrieval.independent_candidate_count,
            "compatible_candidate_count": len(retrieval.pool),
            "independent_compatible_candidate_count": (
                retrieval.independent_compatible_candidate_count
            ),
            "retrieval_trace": tuple(asdict(value) for value in retrieval.trace),
            "recommended_recipe_core_ids": recommended_core_ids,
            "top1_recipe_match": top1_match,
            "topk_recipe_match": topk_match,
            "error": None if recommendations else retrieval.level,
        }
        cases.append(case)
        if recommendations:
            precedent_rows = tuple(value[0] for value in retrieval.pool[:3])
            review_rows.append(
                {
                    "case_id": row.observation_id or row.reaction_id,
                    "query_reaction_smiles": row.reaction_smiles,
                    "query_reaction_label": _reaction_label(row.reaction_label),
                    "query_quality_status": _quality_status(core),
                    "retrieval_level": base_level,
                    "independent_compatible_candidate_count": (
                        retrieval.independent_compatible_candidate_count
                    ),
                    "precedent_reaction_smiles": " || ".join(
                        value.reaction_smiles for value in precedent_rows
                    ),
                    "precedent_reaction_labels": " || ".join(
                        _reaction_label(value.reaction_label)
                        for value in precedent_rows
                    ),
                    "chemist_core_match": "",
                    "chemist_notes": "",
                }
            )

    query_count = len(cases)
    covered = sum(bool(case["valid"]) for case in cases)
    candidate_count = sum(int(case["candidate_count"]) for case in cases)
    compatible_count = sum(
        int(case["compatible_candidate_count"]) for case in cases
    )
    report = {
        "schema_version": REACTION_CORE_EVALUATION_SCHEMA_VERSION,
        "evaluator_version": REACTION_CORE_EVALUATOR_VERSION,
        "records_path": str(Path(records_path)),
        "definition_versions": {
            "reaction_core_retrieval": str(
                load_reaction_core_retrieval_rules()["definition_id"]
            ),
            "reaction_core_quality": str(
                load_reaction_core_quality_rules()["definition_version"]
            ),
            "evidence_support": str(
                load_evidence_support_rules()["definition_version"]
            ),
            "fallback_retrieval": str(
                load_fallback_retrieval_rules()["schema_version"]
            ),
            "ranking": str(load_generic_ranking_rules()["schema_version"]),
        },
        "parameters": {
            "test_fraction": test_fraction,
            "seed": seed,
            "top_k": top_k,
            "minimum_pool_size": minimum_pool_size,
            "split_mode": split_mode,
        },
        "split": {
            "eligible_record_count": len(full_index.rows),
            "train_record_count": len(split.train_rows),
            "test_record_count": len(split.test_rows),
            "train_group_count": len(split.train_group_ids),
            "test_group_count": len(split.test_group_ids),
            "leakage_group_count": len(
                set(split.train_group_ids) & set(split.test_group_ids)
            ),
        },
        "mapping_equivalence": {
            "all_rows": _mapping_summary(full_index.rows),
            "train_rows": _mapping_summary(split.train_rows),
            "test_rows": _mapping_summary(split.test_rows),
        },
        "metrics": {
            "test_record_count": len(split.test_rows),
            "core_query_count": query_count,
            "core_unavailable_count": unavailable,
            "core_availability_rate": round(
                query_count / len(split.test_rows), 6
            )
            if split.test_rows
            else 0.0,
            "covered_query_count": covered,
            "coverage_rate": round(covered / query_count, 6)
            if query_count
            else 0.0,
            "abstention_count": query_count - covered,
            "top1_recipe_recovery_rate": round(
                sum(bool(case["top1_recipe_match"]) for case in cases)
                / query_count,
                6,
            )
            if query_count
            else 0.0,
            "topk_recipe_recovery_rate": round(
                sum(bool(case["topk_recipe_match"]) for case in cases)
                / query_count,
                6,
            )
            if query_count
            else 0.0,
            "candidate_count": candidate_count,
            "compatible_candidate_count": compatible_count,
            "compatibility_survival_rate": round(
                compatible_count / candidate_count,
                6,
            )
            if candidate_count
            else 0.0,
        },
        "by_retrieval_level": _stratified_metrics(cases, "retrieval_level"),
        "by_quality_status": _stratified_metrics(cases, "quality_status"),
        "error_counts": dict(
            sorted(Counter(str(case["error"]) for case in cases if case["error"]).items())
        ),
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "reaction_core_calibration_cases.jsonl").write_text(
        "".join(
            json.dumps(case, ensure_ascii=False, sort_keys=True) + "\n"
            for case in cases
        ),
        encoding="utf-8",
    )
    (destination / "reaction_core_calibration_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (destination / "reaction_core_calibration_report.md").write_text(
        _markdown(report),
        encoding="utf-8",
    )
    with (destination / "reaction_core_chemist_review.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=_REVIEW_FIELDS)
        writer.writeheader()
        writer.writerows(review_rows)
    return report


__all__ = [
    "REACTION_CORE_EVALUATION_SCHEMA_VERSION",
    "REACTION_CORE_EVALUATOR_VERSION",
    "evaluate_reaction_core_index",
]
