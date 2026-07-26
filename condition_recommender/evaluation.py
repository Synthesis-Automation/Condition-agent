"""Leakage-safe offline evaluation for generic reaction recommendation."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

from .compatibility import assess_recipe_compatibility, load_compatibility_rules
from .generic_api import recommend_indexed_signature
from .generic_indexing import (
    GenericIndexedReaction,
    build_generic_index_from_rows,
    load_generic_index,
)
from .generic_retrieval import load_generic_retrieval_rules
from .recipe_ranking import load_generic_ranking_rules
from .similarity import load_generic_similarity_rules


@dataclass(frozen=True)
class GroupedHoldoutSplit:
    """Deterministic split with disjoint canonical-reaction groups."""

    train_rows: Tuple[GenericIndexedReaction, ...]
    test_rows: Tuple[GenericIndexedReaction, ...]
    train_group_ids: Tuple[str, ...]
    test_group_ids: Tuple[str, ...]


def _evaluation_groups(
    rows: Tuple[GenericIndexedReaction, ...],
) -> Dict[str, list[GenericIndexedReaction]]:
    """Join rows connected by either publication or canonical reaction."""
    parents = list(range(len(rows)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    owners: Dict[tuple[str, str], int] = {}
    row_tokens: list[tuple[str, ...]] = []
    for index, row in enumerate(rows):
        tokens = []
        if row.reference_id:
            tokens.append(f"reference:{row.reference_id}")
        if row.canonical_reaction_id:
            tokens.append(f"reaction:{row.canonical_reaction_id}")
        if not tokens:
            tokens.append(
                "observation:"
                + (row.observation_id or row.reaction_id or f"row-{index}")
            )
        row_tokens.append(tuple(tokens))
        for token in tokens:
            key = tuple(token.split(":", maxsplit=1))
            previous = owners.setdefault(key, index)
            union(index, previous)

    component_tokens: Dict[int, set[str]] = {}
    for index, tokens in enumerate(row_tokens):
        component_tokens.setdefault(find(index), set()).update(tokens)
    component_ids = {
        root: "EG1:"
        + hashlib.sha256(
            "\0".join(sorted(tokens)).encode("utf-8")
        ).hexdigest()
        for root, tokens in component_tokens.items()
    }
    groups: Dict[str, list[GenericIndexedReaction]] = {}
    for index, row in enumerate(rows):
        groups.setdefault(component_ids[find(index)], []).append(row)
    return groups


def grouped_holdout_split(
    rows: Iterable[GenericIndexedReaction],
    *,
    test_fraction: float = 0.2,
    seed: int = 17,
) -> GroupedHoldoutSplit:
    """Split whole canonical reactions so duplicates cannot cross the boundary."""
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be between zero and one")
    row_values = tuple(rows)
    groups = _evaluation_groups(row_values)
    if len(groups) < 2:
        raise ValueError("At least two canonical reaction groups are required")

    def order(group_id: str) -> tuple[str, str]:
        digest = hashlib.sha256(f"{seed}\0{group_id}".encode("utf-8")).hexdigest()
        return digest, group_id

    ordered_groups = sorted(groups, key=order)
    test_group_count = min(
        len(ordered_groups) - 1,
        max(1, round(len(ordered_groups) * test_fraction)),
    )
    test_groups = set(ordered_groups[:test_group_count])
    train_rows = []
    test_rows = []
    for group_id in sorted(groups):
        target = test_rows if group_id in test_groups else train_rows
        target.extend(groups[group_id])
    def row_key(row: GenericIndexedReaction) -> tuple[str, str, str, str]:
        return (
            row.reference_id or row.canonical_reaction_id,
            row.reaction_id,
            row.observation_id,
            row.recipe_id,
        )
    return GroupedHoldoutSplit(
        train_rows=tuple(sorted(train_rows, key=row_key)),
        test_rows=tuple(sorted(test_rows, key=row_key)),
        train_group_ids=tuple(sorted(set(groups) - test_groups)),
        test_group_ids=tuple(sorted(test_groups)),
    )


def _mean(values: list[float]) -> float | None:
    return round(sum(values) / len(values), 6) if values else None


def _markdown(report: Dict[str, Any]) -> str:
    metrics = report["metrics"]
    split = report["split"]
    return "\n".join(
        (
            "# Generic recommendation held-out evaluation",
            "",
            f"- Train groups: {split['train_group_count']}",
            f"- Test groups: {split['test_group_count']}",
            f"- Leakage groups: {split['leakage_group_count']}",
            f"- Query coverage: {metrics['coverage_rate']:.2%}",
            f"- Top-1 recipe recovery: {metrics['top1_recipe_recovery_rate']:.2%}",
            f"- Top-k recipe recovery: {metrics['topk_recipe_recovery_rate']:.2%}",
            f"- Yield MAE: {metrics['yield_mae_pct']}",
            f"- Hard-incompatible recommendations: {metrics['hard_incompatible_recommendation_count']}",
            "",
            "Recipe recovery includes all held-out observations. The report also",
            "provides recovery conditional on the recipe being present in training.",
            "",
        )
    )


def evaluate_generic_index(
    records_path: str | Path,
    output_dir: str | Path,
    *,
    test_fraction: float = 0.2,
    seed: int = 17,
    top_k: int = 5,
    minimum_pool_size: int | None = None,
) -> Dict[str, Any]:
    """Evaluate retrieval with canonical-reaction-group holdout protection."""
    if top_k < 1:
        raise ValueError("top_k must be positive")
    full_index = load_generic_index(records_path)
    split = grouped_holdout_split(
        full_index.rows, test_fraction=test_fraction, seed=seed
    )
    train_index = build_generic_index_from_rows(split.train_rows)
    train_recipe_core_ids = {row.recipe_core_id for row in split.train_rows}
    retrieval_levels = Counter()
    errors = Counter()
    cases = []
    covered = top1_matches = topk_matches = 0
    seen_recipe_queries = seen_top1_matches = seen_topk_matches = 0
    excluded_candidates = hard_incompatible = 0
    yield_errors = []
    for row in split.test_rows:
        result = recommend_indexed_signature(
            row.signature,
            train_index,
            query_reaction_smiles=row.reaction_smiles,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
        )
        retrieval_levels[result.retrieval_level or "none"] += 1
        excluded_candidates += result.excluded_candidate_count
        recommended_ids = tuple(item.recipe_id for item in result.recommendations)
        recommended_core_ids = tuple(
            item.recipe_core_id for item in result.recommendations
        )
        recipe_seen = row.recipe_core_id in train_recipe_core_ids
        if recipe_seen:
            seen_recipe_queries += 1
        top1_match = bool(
            recommended_core_ids
            and recommended_core_ids[0] == row.recipe_core_id
        )
        topk_match = row.recipe_core_id in recommended_core_ids
        if result.valid and result.recommendations:
            covered += 1
            top1_matches += int(top1_match)
            topk_matches += int(topk_match)
            seen_top1_matches += int(recipe_seen and top1_match)
            seen_topk_matches += int(recipe_seen and topk_match)
            predicted_yield = result.recommendations[0].expected_yield_pct
            if predicted_yield is not None and row.yield_pct is not None:
                yield_errors.append(abs(predicted_yield - row.yield_pct))
            for recommendation in result.recommendations:
                assessment = assess_recipe_compatibility(
                    row.signature, recommendation.resolved_recipe
                )
                hard_incompatible += int(not assessment.compatible)
        else:
            errors[result.error or "UNKNOWN_ERROR"] += 1
        cases.append(
            {
                "canonical_reaction_id": row.canonical_reaction_id,
                "reaction_id": row.reaction_id,
                "observation_id": row.observation_id,
                "actual_recipe_id": row.recipe_id,
                "actual_recipe_core_id": row.recipe_core_id,
                "actual_yield_pct": row.yield_pct,
                "recipe_seen_in_training": recipe_seen,
                "valid": result.valid,
                "retrieval_level": result.retrieval_level,
                "candidate_count": result.candidate_count,
                "independent_candidate_count": (
                    result.independent_candidate_count
                ),
                "compatible_candidate_count": result.compatible_candidate_count,
                "independent_compatible_candidate_count": (
                    result.independent_compatible_candidate_count
                ),
                "excluded_candidate_count": result.excluded_candidate_count,
                "retrieval_trace": tuple(
                    asdict(item) for item in result.retrieval_trace
                ),
                "recommended_recipe_ids": recommended_ids,
                "recommended_recipe_core_ids": recommended_core_ids,
                "top1_recipe_match": top1_match,
                "topk_recipe_match": topk_match,
                "predicted_yield_pct": result.recommendations[0].expected_yield_pct
                if result.recommendations
                else None,
                "top_recommendation_score_trace": (
                    asdict(result.recommendations[0].score_trace)
                    if result.recommendations
                    else None
                ),
                "error": result.error,
            }
        )
    query_count = len(split.test_rows)
    train_groups = set(split.train_group_ids)
    test_groups = set(split.test_group_ids)
    report: Dict[str, Any] = {
        "schema_version": "1.2",
        "evaluator_version": "generic_grouped_holdout.v1.2",
        "records_path": str(Path(records_path)),
        "definition_versions": {
            "compatibility": str(load_compatibility_rules()["schema_version"]),
            "retrieval": str(load_generic_retrieval_rules()["schema_version"]),
            "similarity": str(load_generic_similarity_rules()["schema_version"]),
            "ranking": str(load_generic_ranking_rules()["schema_version"]),
        },
        "parameters": {
            "test_fraction": test_fraction,
            "seed": seed,
            "top_k": top_k,
            "minimum_pool_size": minimum_pool_size,
        },
        "split": {
            "eligible_record_count": len(full_index.rows),
            "canonical_group_count": len(train_groups | test_groups),
            "train_record_count": len(split.train_rows),
            "test_record_count": len(split.test_rows),
            "train_group_count": len(train_groups),
            "test_group_count": len(test_groups),
            "leakage_group_count": len(train_groups & test_groups),
        },
        "metrics": {
            "query_count": query_count,
            "covered_query_count": covered,
            "coverage_rate": round(covered / query_count, 6) if query_count else 0.0,
            "top1_recipe_match_count": top1_matches,
            "top1_recipe_recovery_rate": round(top1_matches / query_count, 6)
            if query_count
            else 0.0,
            "topk_recipe_match_count": topk_matches,
            "topk_recipe_recovery_rate": round(topk_matches / query_count, 6)
            if query_count
            else 0.0,
            "seen_recipe_query_count": seen_recipe_queries,
            "seen_recipe_top1_recovery_rate": round(
                seen_top1_matches / seen_recipe_queries, 6
            )
            if seen_recipe_queries
            else 0.0,
            "seen_recipe_topk_recovery_rate": round(
                seen_topk_matches / seen_recipe_queries, 6
            )
            if seen_recipe_queries
            else 0.0,
            "yield_mae_pct": _mean(yield_errors),
            "excluded_candidate_count": excluded_candidates,
            "hard_incompatible_recommendation_count": hard_incompatible,
        },
        "retrieval_level_counts": dict(sorted(retrieval_levels.items())),
        "error_counts": dict(sorted(errors.items())),
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "evaluation_cases.jsonl").write_text(
        "".join(
            json.dumps(case, ensure_ascii=False, sort_keys=True) + "\n"
            for case in cases
        ),
        encoding="utf-8",
    )
    (destination / "evaluation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (destination / "evaluation_report.md").write_text(
        _markdown(report), encoding="utf-8"
    )
    return report


__all__ = [
    "GroupedHoldoutSplit",
    "evaluate_generic_index",
    "grouped_holdout_split",
]
