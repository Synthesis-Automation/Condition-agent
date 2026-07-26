"""Leakage-safe offline evaluation for generic reaction recommendation."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Literal, Mapping, Optional, Tuple

from .compatibility import assess_recipe_compatibility, load_compatibility_rules
from .generic_api import recommend_indexed_signature
from .generic_indexing import (
    GenericIndexedReaction,
    build_generic_index_from_rows,
    load_generic_index,
)
from .generic_retrieval import load_generic_retrieval_rules
from .generic_retrieval import RetrievalStrategy
from .recipe_ranking import load_generic_ranking_rules
from .similarity import load_generic_similarity_rules


@dataclass(frozen=True)
class GroupedHoldoutSplit:
    """Deterministic split with explicit leakage-control metadata."""

    train_rows: Tuple[GenericIndexedReaction, ...]
    test_rows: Tuple[GenericIndexedReaction, ...]
    train_group_ids: Tuple[str, ...]
    test_group_ids: Tuple[str, ...]
    split_mode: str = "grouped_random"
    cutoff_year: Optional[int] = None
    omitted_row_count: int = 0


def _evaluation_groups(
    rows: Tuple[GenericIndexedReaction, ...],
    *,
    include_scaffolds: bool = False,
    include_sources: bool = False,
) -> Dict[str, list[GenericIndexedReaction]]:
    """Join rows connected by publication, reaction, and optional scaffolds."""
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
        if include_scaffolds:
            tokens.extend(
                f"scaffold:{token}" for token in row.scaffold_tokens
            )
        if include_sources and row.source_dataset:
            tokens.append(f"source:{row.source_dataset}")
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
        root: (
            "ESG1:"
            if include_scaffolds
            else "ECG1:"
            if include_sources
            else "EG1:"
        )
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
    split_mode: Literal[
        "grouped_random",
        "scaffold_disjoint",
        "source_disjoint",
        "forward_time",
    ] = "grouped_random",
) -> GroupedHoldoutSplit:
    """Split connected evidence groups under the requested leakage policy."""
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be between zero and one")
    row_values = tuple(rows)
    if split_mode == "forward_time":
        return _forward_time_split(row_values, test_fraction=test_fraction)
    if split_mode not in {
        "grouped_random",
        "scaffold_disjoint",
        "source_disjoint",
    }:
        raise ValueError(f"Unsupported evaluation split mode: {split_mode}")
    groups = _evaluation_groups(
        row_values,
        include_scaffolds=split_mode == "scaffold_disjoint",
        include_sources=split_mode == "source_disjoint",
    )
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
        split_mode=split_mode,
    )


def _forward_time_split(
    rows: Tuple[GenericIndexedReaction, ...],
    *,
    test_fraction: float,
) -> GroupedHoldoutSplit:
    """Hold out the latest connected groups without future-to-past leakage."""
    groups = _evaluation_groups(rows)
    dated_groups = []
    omitted = 0
    for group_id, members in groups.items():
        years = tuple(
            row.publication_year
            for row in members
            if row.publication_year is not None
        )
        if not years:
            omitted += len(members)
            continue
        dated_groups.append((max(years), group_id, members))
    if len(dated_groups) < 2:
        raise ValueError("At least two dated evaluation groups are required")
    dated_groups.sort(key=lambda item: (item[0], item[1]))
    test_group_count = min(
        len(dated_groups) - 1,
        max(1, round(len(dated_groups) * test_fraction)),
    )
    cutoff_year = dated_groups[-test_group_count][0]
    test_values = [
        item for item in dated_groups if item[0] >= cutoff_year
    ]
    train_values = [
        item for item in dated_groups if item[0] < cutoff_year
    ]
    if not train_values or not test_values:
        raise ValueError("Forward-time split requires at least two distinct years")

    def row_key(row: GenericIndexedReaction) -> tuple[str, str, str, str]:
        return (
            row.reference_id or row.canonical_reaction_id,
            row.reaction_id,
            row.observation_id,
            row.recipe_id,
        )

    return GroupedHoldoutSplit(
        train_rows=tuple(
            sorted(
                (
                    row
                    for _, _, members in train_values
                    for row in members
                ),
                key=row_key,
            )
        ),
        test_rows=tuple(
            sorted(
                (
                    row
                    for _, _, members in test_values
                    for row in members
                ),
                key=row_key,
            )
        ),
        train_group_ids=tuple(sorted(item[1] for item in train_values)),
        test_group_ids=tuple(sorted(item[1] for item in test_values)),
        split_mode="forward_time",
        cutoff_year=cutoff_year,
        omitted_row_count=omitted,
    )


def _leakage_summary(split: GroupedHoldoutSplit) -> Dict[str, Any]:
    train_references = {row.reference_id for row in split.train_rows if row.reference_id}
    test_references = {row.reference_id for row in split.test_rows if row.reference_id}
    train_reactions = {
        row.canonical_reaction_id
        for row in split.train_rows
        if row.canonical_reaction_id
    }
    test_reactions = {
        row.canonical_reaction_id
        for row in split.test_rows
        if row.canonical_reaction_id
    }
    train_scaffolds = {
        token for row in split.train_rows for token in row.scaffold_tokens
    }
    test_scaffolds = {
        token for row in split.test_rows for token in row.scaffold_tokens
    }
    train_sources = {row.source_dataset for row in split.train_rows}
    test_sources = {row.source_dataset for row in split.test_rows}
    train_years = sorted(
        row.publication_year
        for row in split.train_rows
        if row.publication_year is not None
    )
    test_years = sorted(
        row.publication_year
        for row in split.test_rows
        if row.publication_year is not None
    )
    return {
        "reference_overlap_count": len(train_references & test_references),
        "canonical_reaction_overlap_count": len(train_reactions & test_reactions),
        "scaffold_overlap_count": len(train_scaffolds & test_scaffolds),
        "source_dataset_overlap_count": len(train_sources & test_sources),
        "train_scaffold_count": len(train_scaffolds),
        "test_scaffold_count": len(test_scaffolds),
        "train_year_range": (
            [train_years[0], train_years[-1]] if train_years else None
        ),
        "test_year_range": (
            [test_years[0], test_years[-1]] if test_years else None
        ),
    }


def _mean(values: list[float]) -> float | None:
    return round(sum(values) / len(values), 6) if values else None


def _stratified_metrics(
    cases: Iterable[Mapping[str, Any]],
    field: str,
) -> Dict[str, Dict[str, Any]]:
    groups: Dict[str, list[Mapping[str, Any]]] = {}
    for case in cases:
        groups.setdefault(str(case.get(field) or "unknown"), []).append(case)
    values = {}
    for key, members in sorted(groups.items()):
        count = len(members)
        covered = sum(
            int(bool(member["valid"] and member["recommended_recipe_core_ids"]))
            for member in members
        )
        seen = [member for member in members if member["recipe_seen_in_training"]]
        values[key] = {
            "query_count": count,
            "coverage_rate": round(covered / count, 6) if count else 0.0,
            "top1_recipe_recovery_rate": round(
                sum(int(member["top1_recipe_match"]) for member in members)
                / count,
                6,
            )
            if count
            else 0.0,
            "topk_recipe_recovery_rate": round(
                sum(int(member["topk_recipe_match"]) for member in members)
                / count,
                6,
            )
            if count
            else 0.0,
            "seen_recipe_query_count": len(seen),
            "seen_top1_recovery_rate": round(
                sum(int(member["top1_recipe_match"]) for member in seen)
                / len(seen),
                6,
            )
            if seen
            else 0.0,
            "hard_incompatible_recommendation_count": sum(
                int(member["hard_incompatible_recommendation_count"])
                for member in members
            ),
        }
    return values


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
    split_mode: Literal[
        "grouped_random",
        "scaffold_disjoint",
        "source_disjoint",
        "forward_time",
    ] = "grouped_random",
    retrieval_strategy: RetrievalStrategy = "hybrid",
    ranking_weights: Mapping[str, float] | None = None,
) -> Dict[str, Any]:
    """Evaluate retrieval with canonical-reaction-group holdout protection."""
    if top_k < 1:
        raise ValueError("top_k must be positive")
    full_index = load_generic_index(records_path)
    split = grouped_holdout_split(
        full_index.rows,
        test_fraction=test_fraction,
        seed=seed,
        split_mode=split_mode,
    )
    train_index = build_generic_index_from_rows(split.train_rows)
    train_recipe_core_ids = {row.recipe_core_id for row in split.train_rows}
    train_recipe_ids = {row.recipe_id for row in split.train_rows}
    retrieval_levels = Counter()
    errors = Counter()
    cases = []
    covered = top1_matches = topk_matches = 0
    seen_recipe_queries = seen_top1_matches = seen_topk_matches = 0
    seen_variant_queries = seen_variant_top1 = seen_variant_topk = 0
    excluded_candidates = hard_incompatible = 0
    explanation_complete = uncertain_recommendations = 0
    independent_support_values = []
    yield_errors = []
    for row in split.test_rows:
        result = recommend_indexed_signature(
            row.signature,
            train_index,
            query_reaction_smiles=row.reaction_smiles,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            retrieval_strategy=retrieval_strategy,
            ranking_weights=ranking_weights,
        )
        retrieval_levels[result.retrieval_level or "none"] += 1
        excluded_candidates += result.excluded_candidate_count
        recommended_ids = tuple(item.recipe_id for item in result.recommendations)
        recommended_core_ids = tuple(
            item.recipe_core_id for item in result.recommendations
        )
        recipe_seen = row.recipe_core_id in train_recipe_core_ids
        variant_seen = row.recipe_id in train_recipe_ids
        if recipe_seen:
            seen_recipe_queries += 1
        if variant_seen:
            seen_variant_queries += 1
        top1_match = bool(
            recommended_core_ids
            and recommended_core_ids[0] == row.recipe_core_id
        )
        topk_match = row.recipe_core_id in recommended_core_ids
        top1_variant_match = bool(
            recommended_ids and recommended_ids[0] == row.recipe_id
        )
        topk_variant_match = row.recipe_id in recommended_ids
        case_hard_incompatible = 0
        if result.valid and result.recommendations:
            covered += 1
            top1_matches += int(top1_match)
            topk_matches += int(topk_match)
            seen_top1_matches += int(recipe_seen and top1_match)
            seen_topk_matches += int(recipe_seen and topk_match)
            seen_variant_top1 += int(variant_seen and top1_variant_match)
            seen_variant_topk += int(variant_seen and topk_variant_match)
            top_recommendation = result.recommendations[0]
            independent_support_values.append(
                top_recommendation.score_trace.independent_evidence_count
            )
            uncertain_recommendations += int(
                (
                    top_recommendation.score_trace.ranking_components.get(
                        "condition_certainty"
                    )
                    or 0.0
                )
                < 1.0
            )
            explanation_complete += int(
                bool(
                    top_recommendation.explanation
                    and top_recommendation.precedent_reaction_ids
                    and top_recommendation.score_trace.definition_versions
                )
            )
            predicted_yield = result.recommendations[0].expected_yield_pct
            if predicted_yield is not None and row.yield_pct is not None:
                yield_errors.append(abs(predicted_yield - row.yield_pct))
            for recommendation in result.recommendations:
                assessment = assess_recipe_compatibility(
                    row.signature, recommendation.resolved_recipe
                )
                incompatible = int(not assessment.compatible)
                hard_incompatible += incompatible
                case_hard_incompatible += incompatible
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
                "transformation_class": row.transformation_class,
                "evidence_quality": str(
                    row.signature.get("evidence_quality")
                    or row.chemistry_status
                    or "unknown"
                ),
                "scaffold_key": row.scaffold_key,
                "recipe_seen_in_training": recipe_seen,
                "variant_seen_in_training": variant_seen,
                "valid": result.valid,
                "retrieval_level": result.retrieval_level,
                "retrieval_strategy": result.retrieval_strategy,
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
                "top1_variant_match": top1_variant_match,
                "topk_variant_match": topk_variant_match,
                "hard_incompatible_recommendation_count": (
                    case_hard_incompatible
                ),
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
        "schema_version": "1.4",
        "evaluator_version": "generic_leakage_safe.v1.4",
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
            "split_mode": split_mode,
            "retrieval_strategy": retrieval_strategy,
            "ranking_weights": (
                dict(sorted(ranking_weights.items()))
                if ranking_weights is not None
                else None
            ),
        },
        "split": {
            "eligible_record_count": len(full_index.rows),
            "canonical_group_count": len(train_groups | test_groups),
            "train_record_count": len(split.train_rows),
            "test_record_count": len(split.test_rows),
            "train_group_count": len(train_groups),
            "test_group_count": len(test_groups),
            "leakage_group_count": len(train_groups & test_groups),
            "split_mode": split.split_mode,
            "cutoff_year": split.cutoff_year,
            "omitted_row_count": split.omitted_row_count,
            **_leakage_summary(split),
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
            "seen_variant_query_count": seen_variant_queries,
            "seen_variant_top1_recovery_rate": round(
                seen_variant_top1 / seen_variant_queries, 6
            )
            if seen_variant_queries
            else 0.0,
            "seen_variant_topk_recovery_rate": round(
                seen_variant_topk / seen_variant_queries, 6
            )
            if seen_variant_queries
            else 0.0,
            "abstention_count": query_count - covered,
            "abstention_rate": round((query_count - covered) / query_count, 6)
            if query_count
            else 0.0,
            "yield_mae_pct": _mean(yield_errors),
            "yield_prediction_count": len(yield_errors),
            "excluded_candidate_count": excluded_candidates,
            "hard_incompatible_recommendation_count": hard_incompatible,
            "condition_uncertain_recommendation_count": (
                uncertain_recommendations
            ),
            "condition_uncertain_recommendation_rate": round(
                uncertain_recommendations / covered, 6
            )
            if covered
            else 0.0,
            "explanation_complete_count": explanation_complete,
            "explanation_complete_rate": round(
                explanation_complete / covered, 6
            )
            if covered
            else 0.0,
            "mean_top_independent_support": _mean(
                [float(value) for value in independent_support_values]
            ),
        },
        "by_transformation_class": _stratified_metrics(
            cases, "transformation_class"
        ),
        "by_evidence_quality": _stratified_metrics(cases, "evidence_quality"),
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


def compare_generic_baselines(
    records_path: str | Path,
    output_dir: str | Path,
    *,
    test_fraction: float = 0.2,
    seed: int = 17,
    top_k: int = 5,
    minimum_pool_size: int | None = None,
    split_mode: Literal[
        "grouped_random",
        "scaffold_disjoint",
        "source_disjoint",
        "forward_time",
    ] = "grouped_random",
) -> Dict[str, Any]:
    """Compare retrieval baselines on one deterministic leakage-safe split."""
    strategies: Tuple[RetrievalStrategy, ...] = (
        "family_only",
        "generic_only",
        "hybrid",
        "transformation_prior",
        "legacy_pilot",
    )
    destination = Path(output_dir)
    reports = {}
    for strategy in strategies:
        reports[strategy] = evaluate_generic_index(
            records_path,
            destination / strategy,
            test_fraction=test_fraction,
            seed=seed,
            top_k=top_k,
            minimum_pool_size=minimum_pool_size,
            split_mode=split_mode,
            retrieval_strategy=strategy,
        )
    metric_names = (
        "coverage_rate",
        "top1_recipe_recovery_rate",
        "topk_recipe_recovery_rate",
        "seen_recipe_top1_recovery_rate",
        "seen_recipe_topk_recovery_rate",
        "yield_mae_pct",
        "hard_incompatible_recommendation_count",
    )
    comparison: Dict[str, Any] = {
        "schema_version": "1.0",
        "artifact_type": "generic_baseline_comparison",
        "records_path": str(Path(records_path)),
        "parameters": {
            "test_fraction": test_fraction,
            "seed": seed,
            "top_k": top_k,
            "minimum_pool_size": minimum_pool_size,
            "split_mode": split_mode,
        },
        "strategies": {
            strategy: {
                name: reports[strategy]["metrics"][name]
                for name in metric_names
            }
            for strategy in strategies
        },
        "split": reports["hybrid"]["split"],
    }
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "baseline_comparison.json").write_text(
        json.dumps(comparison, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    markdown = [
        "# Generic recommendation baseline comparison",
        "",
        "| Strategy | Coverage | Top-1 | Top-k | Seen top-1 | Hard conflicts |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for strategy in strategies:
        metrics = comparison["strategies"][strategy]
        markdown.append(
            f"| {strategy} | {metrics['coverage_rate']:.2%} | "
            f"{metrics['top1_recipe_recovery_rate']:.2%} | "
            f"{metrics['topk_recipe_recovery_rate']:.2%} | "
            f"{metrics['seen_recipe_top1_recovery_rate']:.2%} | "
            f"{metrics['hard_incompatible_recommendation_count']} |"
        )
    (destination / "baseline_comparison.md").write_text(
        "\n".join(markdown) + "\n",
        encoding="utf-8",
    )
    return comparison


__all__ = [
    "GroupedHoldoutSplit",
    "compare_generic_baselines",
    "evaluate_generic_index",
    "grouped_holdout_split",
]
