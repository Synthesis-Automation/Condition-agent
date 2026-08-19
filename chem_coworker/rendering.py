"""Deterministic text rendering for typed recommendation results."""

from __future__ import annotations

import json
from typing import Any, Iterable, Mapping

from condition_recommender import GenericRecommendationResult

from .contracts import ConditionReview


def _reaction_label(result: GenericRecommendationResult) -> str:
    label = result.reaction_label or {}
    return str(label.get("concise") or label.get("detailed") or "Unlabeled reaction")


def _component_text(recipe: Mapping[str, Any]) -> str:
    components = recipe.get("components")
    if not isinstance(components, Iterable) or isinstance(components, (str, bytes)):
        return str(
            recipe.get("recipe_core_id") or recipe.get("recipe_id") or "resolved recipe"
        )
    rendered = []
    for component in components:
        if not isinstance(component, Mapping):
            continue
        name = (
            component.get("canonical_name")
            or component.get("raw_identifier")
            or component.get("substance_id")
        )
        if not name:
            continue
        role = component.get("primary_role")
        rendered.append(f"{name} ({role})" if role else str(name))
    if rendered:
        return ", ".join(rendered)
    return json.dumps(dict(recipe), ensure_ascii=False, sort_keys=True)


def _ordered_recommendations(
    result: GenericRecommendationResult,
    review: ConditionReview | None,
) -> tuple[Any, ...]:
    if review is None or not review.presentation_recipe_ids:
        return tuple(result.recommendations)
    by_id = {item.recipe_id: item for item in result.recommendations}
    ordered = [
        by_id[recipe_id]
        for recipe_id in review.presentation_recipe_ids
        if recipe_id in by_id
    ]
    included = {item.recipe_id for item in ordered}
    ordered.extend(
        item for item in result.recommendations if item.recipe_id not in included
    )
    return tuple(ordered)


def render_recommendation(
    result: GenericRecommendationResult,
    review: ConditionReview | None = None,
) -> str:
    """Render a concise answer while retaining the typed result as source of truth."""

    if not result.valid:
        lines = [f"No condition recommendation: {result.error or 'invalid reaction'}."]
        proposal = result.completion_proposal or {}
        if proposal.get("requirements"):
            lines.append(
                "Confirm the proposed missing reaction source(s) and resubmit."
            )
        if result.warnings:
            lines.append("Warnings: " + "; ".join(result.warnings))
        return "\n".join(lines)

    heading = _reaction_label(result)
    family = result.named_family or "unassigned (generic structural retrieval)"
    lines = [
        f"Reaction: {heading}",
        f"Interpretation: {result.transformation_class or 'unresolved'}; family: {family}",
        f"Retrieval: {result.retrieval_level or 'none'} ({result.candidate_count} candidates)",
    ]
    if not result.recommendations:
        lines.append("No compatible condition recipe met the recommendation criteria.")
    else:
        lines.append("Recommended conditions:")
        reviews_by_id = {
            item.recipe_id: item for item in (review.candidates if review else ())
        }
        for display_rank, item in enumerate(
            _ordered_recommendations(result, review), start=1
        ):
            yield_text = (
                f", expected yield {item.expected_yield_pct:.1f}%"
                if item.expected_yield_pct is not None
                else ""
            )
            lines.append(
                f"{display_rank}. {_component_text(item.resolved_recipe)} "
                f"[score {item.score:.3f}, {item.reference_support} references{yield_text}]"
            )
            candidate_review = reviews_by_id.get(item.recipe_id)
            if candidate_review is not None:
                lines.append(
                    "   LLM review: "
                    f"{candidate_review.verdict} "
                    f"({candidate_review.confidence:.0%}) — "
                    f"{candidate_review.rationale}"
                )
            if item.explanation:
                lines.append("   Evidence: " + "; ".join(item.explanation))
            if item.cautions:
                lines.append("   Cautions: " + "; ".join(item.cautions))
            if item.precedent_reaction_ids:
                lines.append("   Precedents: " + ", ".join(item.precedent_reaction_ids))
    if result.warnings:
        lines.append("Warnings: " + "; ".join(result.warnings))
    if review is not None:
        if review.status == "completed":
            lines.append(f"LLM review: {review.summary} [{review.model}; advisory]")
            if review.questions:
                lines.append("Questions: " + "; ".join(review.questions))
        elif review.status == "failed" and review.warning:
            lines.append("LLM review unavailable: " + review.warning)
    return "\n".join(lines)
