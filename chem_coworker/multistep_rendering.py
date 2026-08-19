"""Deterministic text rendering for multistep route-search responses."""

from __future__ import annotations

from typing import Sequence

from core_retrosynthesis import MultistepRetrosynthesisRoute

from .contracts import MultistepRetrosynthesisResponse


def ordered_multistep_routes(
    response: MultistepRetrosynthesisResponse,
) -> tuple[MultistepRetrosynthesisRoute, ...]:
    """Return the reviewed route set in advisory presentation order."""

    if response.result is None:
        return ()
    source: Sequence[MultistepRetrosynthesisRoute]
    source = response.result.routes or response.result.partial_routes
    review = response.review
    if review is None or not review.presentation_route_ids:
        return tuple(source)
    by_id = {route.route_id: route for route in source}
    ordered = [
        by_id[route_id]
        for route_id in review.presentation_route_ids
        if route_id in by_id
    ]
    included = {route.route_id for route in ordered}
    ordered.extend(route for route in source if route.route_id not in included)
    return tuple(ordered)


def render_multistep_retrosynthesis(
    response: MultistepRetrosynthesisResponse,
) -> str:
    """Render a concise route summary while retaining typed evidence as truth."""

    if not response.valid or response.result is None:
        return "No multistep result: " + (response.error or "invalid target")
    result = response.result
    route_kind = "solved" if result.routes else "partial"
    lines = [
        f"Target: {response.request.target_smiles}",
        (
            f"Routes: {len(result.routes)} solved, "
            f"{len(result.partial_routes)} partial; showing {route_kind} routes"
        ),
        (
            f"Search: depth <= {result.max_depth}; "
            f"{result.diagnostics.expanded_states} states expanded; "
            f"{result.diagnostics.one_step_calls} one-step calls"
        ),
    ]
    reviews = {
        item.route_id: item
        for item in (response.review.routes if response.review else ())
    }
    for shown_rank, route in enumerate(
        ordered_multistep_routes(response)[: response.request.top_k], start=1
    ):
        original_source = result.routes or result.partial_routes
        original_rank = next(
            index
            for index, value in enumerate(original_source, start=1)
            if value.route_id == route.route_id
        )
        lines.append(
            f"{shown_rank}. {route_kind.title()} route (original rank "
            f"{original_rank}): {route.reaction_count} steps, cost "
            f"{route.route_cost:.3f}"
        )
        for step_number, step in enumerate(route.steps, start=1):
            condition = step.condition_evidence
            condition_text = condition.status if condition else "not evaluated"
            lines.append(
                f"   Step {step_number}: {'.'.join(step.precursor_smiles)} "
                f">> {step.product_smiles} [{condition_text}]"
            )
        unresolved = [leaf.canonical_smiles for leaf in route.leaves if not leaf.terminal]
        if unresolved:
            lines.append("   Unresolved leaves: " + ", ".join(unresolved))
        review = reviews.get(route.route_id)
        if review is not None:
            lines.append(
                f"   LLM review: {review.verdict} ({review.confidence:.0%}) — "
                f"{review.rationale}"
            )
    if response.review is not None:
        if response.review.status == "completed":
            lines.append(
                f"LLM route review: {response.review.summary} "
                f"[{response.review.model}; advisory]"
            )
            if response.review.questions:
                lines.append("Questions: " + "; ".join(response.review.questions))
        elif response.review.status == "failed" and response.review.warning:
            lines.append("LLM review unavailable: " + response.review.warning)
    if response.warnings:
        lines.append("Warnings: " + "; ".join(response.warnings))
    return "\n".join(lines)


__all__ = ["ordered_multistep_routes", "render_multistep_retrosynthesis"]
