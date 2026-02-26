"""Source planning for condition recommendation requests."""

from __future__ import annotations

from typing import List

from .models import (
    OutputView,
    RecommendationRequest,
    RunStrategy,
    SourceGroup,
    SourcePlan,
)


_CORE_SOURCE_ORDER = (
    SourceGroup.LITERATURE.value,
    SourceGroup.EXPERIMENTS.value,
    SourceGroup.RULES.value,
)


def plan_sources(request: RecommendationRequest) -> SourcePlan:
    """Plan which sources to execute and whether dataset loading is required."""
    strategy = request.normalized_run_strategy()
    output_view = request.normalized_output_view()
    source = request.normalized_source_group()

    notes: List[str] = []
    if strategy is RunStrategy.ANALYSIS_ONLY:
        return SourcePlan(
            run_strategy=strategy,
            output_view=output_view,
            notes=("analysis_only; no dataset-backed recommendation execution",),
        )

    needs_precedent = output_view is OutputView.PRECEDENT_ONLY
    if output_view is OutputView.PRECEDENT_ONLY and source in {SourceGroup.ANY, SourceGroup.EXPERIMENTS, SourceGroup.RULES}:
        notes.append("precedent_only output requires literature source execution")

    if strategy is RunStrategy.PER_SOURCE:
        if source is SourceGroup.ANY:
            sources = list(_CORE_SOURCE_ORDER)
        else:
            sources = [source.value]
        if needs_precedent and SourceGroup.LITERATURE.value not in sources:
            sources.insert(0, SourceGroup.LITERATURE.value)
        return SourcePlan(
            sources_to_run=tuple(sources),
            needs_hte_data=bool(sources),
            needs_precedent_data=needs_precedent,
            run_strategy=strategy,
            output_view=output_view,
            notes=tuple(notes),
        )

    # Single-pass
    single_source = None if source is SourceGroup.ANY else source.value
    if needs_precedent and single_source != SourceGroup.LITERATURE.value:
        single_source = SourceGroup.LITERATURE.value
        notes.append("single-pass precedent_only coerced source_group to literature")
    return SourcePlan(
        sources_to_run=(() if single_source is None else (single_source,)),
        single_run_source_group=single_source,
        needs_hte_data=True,
        needs_precedent_data=needs_precedent,
        run_strategy=strategy,
        output_view=output_view,
        notes=tuple(notes),
    )
