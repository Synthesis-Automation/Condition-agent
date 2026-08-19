"""Deterministic text rendering for one-step retrosynthesis responses."""

from __future__ import annotations

from typing import Any, Mapping

from .contracts import RetrosynthesisResponse
from .rendering import _component_text


def ordered_retrosynthesis_strategies(response: RetrosynthesisResponse):
    """Return advisory presentation order without mutating source ranking."""

    review = response.review
    if review is None or not review.presentation_strategy_ids:
        return response.strategies
    by_id = {item.strategy_id: item for item in response.strategies}
    ordered = [
        by_id[strategy_id]
        for strategy_id in review.presentation_strategy_ids
        if strategy_id in by_id
    ]
    included = {item.strategy_id for item in ordered}
    ordered.extend(
        item for item in response.strategies if item.strategy_id not in included
    )
    return tuple(ordered)


def _best_condition(evidence: Any) -> str:
    if evidence is None:
        return "not evaluated"
    if not evidence.recommendations:
        return evidence.status.replace("_", " ")
    recipe = evidence.recommendations[0]
    resolved = recipe.get("resolved_recipe")
    if isinstance(resolved, Mapping):
        return _component_text(resolved)
    return str(recipe.get("recipe_core_id") or evidence.status).replace("_", " ")


def render_retrosynthesis(response: RetrosynthesisResponse) -> str:
    """Render concise strategies while retaining the typed result as truth."""

    if not response.valid:
        return "No retrosynthesis result: " + (response.error or "invalid target")
    lines = [
        f"Target: {response.request.target_smiles}",
        (
            f"Validated strategies: {len(response.strategies)} in review pool; "
            f"showing up to {response.request.top_k}"
        ),
    ]
    if not response.strategies:
        lines.append("No forward-validated one-step disconnection was found.")
    else:
        lines.append("One-step retrosynthesis strategies:")
        review_by_id = {
            item.strategy_id: item
            for item in (response.review.candidates if response.review else ())
        }
        conditions_by_id = {
            item.strategy_id: item.evidence for item in response.condition_evidence
        }
        ordered = ordered_retrosynthesis_strategies(response)[
            : response.request.top_k
        ]
        for display_rank, strategy in enumerate(ordered, start=1):
            candidate = strategy.representative
            lines.append(
                f"{display_rank}. {candidate.precursor_smiles} "
                f"[score {candidate.score:.3f}, {candidate.abstraction_level}, "
                f"{strategy.independent_reference_support} references]"
            )
            lines.append(
                "   Forward check: " + candidate.forward_validation_status
            )
            if candidate.transformation_kind:
                lines.append(
                    "   Transformation: " + candidate.transformation_kind
                )
            if strategy.total_realization_count > 1:
                lines.append(
                    "   Concrete variants: "
                    f"{strategy.total_realization_count} total; "
                    f"{len(strategy.realizations)} retained"
                )
            condition = conditions_by_id.get(strategy.strategy_id)
            if response.request.include_conditions:
                lines.append("   Conditions: " + _best_condition(condition))
            model_review = review_by_id.get(strategy.strategy_id)
            if model_review is not None:
                lines.append(
                    "   LLM review: "
                    f"{model_review.verdict} ({model_review.confidence:.0%}) — "
                    f"{model_review.rationale}"
                )
            if candidate.selectivity_warnings:
                lines.append(
                    "   Selectivity cautions: "
                    + "; ".join(
                        str(item.message) for item in candidate.selectivity_warnings
                    )
                )
            if candidate.precursor_compatibility_disposition != "pass":
                lines.append(
                    "   Precursor compatibility: "
                    + candidate.precursor_compatibility_disposition
                )
            if strategy.precedent_reaction_ids:
                lines.append(
                    "   Precedents: "
                    + ", ".join(strategy.precedent_reaction_ids[:5])
                )
    if response.review is not None:
        if response.review.status == "completed":
            lines.append(
                f"LLM review: {response.review.summary} "
                f"[{response.review.model}; advisory]"
            )
            if response.review.questions:
                lines.append("Questions: " + "; ".join(response.review.questions))
        elif response.review.status == "failed" and response.review.warning:
            lines.append("LLM review unavailable: " + response.review.warning)
    if response.warnings:
        lines.append("Warnings: " + "; ".join(response.warnings))
    return "\n".join(lines)


__all__ = ["ordered_retrosynthesis_strategies", "render_retrosynthesis"]
