"""Deterministic application of registry-normalized condition constraints."""

from __future__ import annotations

from dataclasses import replace
from typing import Dict, List

from condition_registry import ConditionConstraintSet, condition_constraint_conflicts

from .models import GenericRecommendationResult


def apply_condition_constraints(
    result: GenericRecommendationResult,
    constraints: ConditionConstraintSet,
    *,
    top_k: int,
) -> GenericRecommendationResult:
    """Remove conflicting recipes before returning the canonical ranked result."""

    if not constraints.constraints:
        return replace(
            result,
            recommendations=result.recommendations[:top_k],
            condition_constraints=constraints.to_dict(),
        )
    accepted = []
    excluded: List[Dict[str, object]] = []
    for recommendation in result.recommendations:
        recipe = dict(recommendation.resolved_recipe)
        if recommendation.synthesis_protocol:
            recipe.setdefault("synthesis_protocol", recommendation.synthesis_protocol)
        conflicts = condition_constraint_conflicts(recipe, constraints)
        if conflicts:
            excluded.append(
                {
                    "recipe_id": recommendation.recipe_id,
                    "constraint_conflicts": list(conflicts),
                }
            )
            continue
        accepted.append(recommendation)
    reranked = tuple(
        replace(item, rank=index)
        for index, item in enumerate(accepted[:top_k], start=1)
    )
    warnings = list(result.warnings)
    warnings.append("CONFIRMED_CONDITION_CONSTRAINTS_APPLIED")
    if excluded:
        warnings.append(f"CONDITION_CONSTRAINT_RECIPES_EXCLUDED:{len(excluded)}")
    valid = result.valid
    error = result.error
    if result.valid and not reranked:
        valid = False
        error = "NO_RECOMMENDATIONS_MATCH_CONFIRMED_CONSTRAINTS"
        warnings.append(error)
    return replace(
        result,
        valid=valid,
        error=error,
        recommendations=reranked,
        condition_constraints=constraints.to_dict(),
        condition_constraint_trace=tuple(excluded),
        warnings=tuple(dict.fromkeys(warnings)),
    )
