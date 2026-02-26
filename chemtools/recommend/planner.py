"""Source planning for condition recommendation requests."""

from __future__ import annotations

from typing import List

from .models import (
    OutputView,
    RecommendationStrategy,
    RecommendationRequest,
    RunStrategy,
    SourceGroup,
    SourcePlan,
)


_CORE_SOURCE_ORDER = (
    SourceGroup.LITERATURE.value,
    SourceGroup.MOTIF.value,
    SourceGroup.RULES.value,
)


def _strategy_plan(
    *,
    run_strategy: RunStrategy,
    output_view: OutputView,
    source: SourceGroup,
    rec_strategy: RecommendationStrategy,
) -> SourcePlan:
    notes: List[str] = []
    if run_strategy is RunStrategy.ANALYSIS_ONLY:
        return SourcePlan(
            run_strategy=run_strategy,
            output_view=output_view,
            notes=("analysis_only; no dataset-backed recommendation execution",),
            recommendation_strategy=rec_strategy.value,
        )

    if rec_strategy is RecommendationStrategy.SIMILARITY:
        if source not in {SourceGroup.ANY, SourceGroup.LITERATURE}:
            notes.append("similarity strategy coerced source_group to literature")
        return SourcePlan(
            sources_to_run=(SourceGroup.LITERATURE.value,),
            single_run_source_group=SourceGroup.LITERATURE.value,
            needs_hte_data=True,
            needs_precedent_data=True,
            run_strategy=RunStrategy.SINGLE_PASS,
            output_view=output_view,
            notes=tuple(notes),
            recommendation_strategy=rec_strategy.value,
        )

    if rec_strategy is RecommendationStrategy.LITERATURE:
        if source not in {SourceGroup.ANY, SourceGroup.LITERATURE}:
            notes.append("literature strategy coerced source_group to literature")
        return SourcePlan(
            sources_to_run=(SourceGroup.LITERATURE.value,),
            single_run_source_group=SourceGroup.LITERATURE.value,
            needs_hte_data=True,
            needs_precedent_data=(output_view is OutputView.PRECEDENT_ONLY),
            run_strategy=RunStrategy.SINGLE_PASS,
            output_view=output_view,
            notes=tuple(notes),
            recommendation_strategy=rec_strategy.value,
        )

    if rec_strategy is RecommendationStrategy.RULES:
        if source not in {SourceGroup.ANY, SourceGroup.RULES}:
            notes.append("rules strategy coerced source_group to rules")
        return SourcePlan(
            sources_to_run=(SourceGroup.RULES.value,),
            single_run_source_group=SourceGroup.RULES.value,
            needs_hte_data=True,
            needs_precedent_data=False,
            run_strategy=RunStrategy.SINGLE_PASS,
            output_view=output_view,
            notes=tuple(notes),
            recommendation_strategy=rec_strategy.value,
        )

    # Motif strategy (motif-based only): experiments + rules by default.
    if source in {SourceGroup.MOTIF, SourceGroup.RULES}:
        planned_run_strategy = (
            RunStrategy.PER_SOURCE if run_strategy is RunStrategy.PER_SOURCE else RunStrategy.SINGLE_PASS
        )
        if planned_run_strategy is RunStrategy.PER_SOURCE:
            sources = (source.value,)
            single_source = None
        else:
            sources = (source.value,)
            single_source = source.value
        return SourcePlan(
            sources_to_run=sources,
            single_run_source_group=single_source,
            needs_hte_data=True,
            needs_precedent_data=False,
            run_strategy=planned_run_strategy,
            output_view=output_view,
            notes=tuple(notes),
            recommendation_strategy=rec_strategy.value,
        )

    if source is SourceGroup.LITERATURE:
        notes.append("motif strategy ignores literature source override")

    return SourcePlan(
        sources_to_run=(SourceGroup.MOTIF.value,),
        single_run_source_group=SourceGroup.MOTIF.value,
        needs_hte_data=True,
        needs_precedent_data=False,
        run_strategy=RunStrategy.SINGLE_PASS,
        output_view=output_view,
        notes=tuple(notes),
        recommendation_strategy=rec_strategy.value,
    )


def plan_sources(request: RecommendationRequest) -> SourcePlan:
    """Plan which sources to execute and whether dataset loading is required."""
    run_strategy = request.normalized_run_strategy()
    output_view = request.normalized_output_view()
    source = request.normalized_source_group()
    rec_strategy = request.normalized_strategy()

    if rec_strategy is not None:
        return _strategy_plan(
            run_strategy=run_strategy,
            output_view=output_view,
            source=source,
            rec_strategy=rec_strategy,
        )

    notes: List[str] = []
    if run_strategy is RunStrategy.ANALYSIS_ONLY:
        return SourcePlan(
            run_strategy=run_strategy,
            output_view=output_view,
            notes=("analysis_only; no dataset-backed recommendation execution",),
        )

    needs_precedent = output_view is OutputView.PRECEDENT_ONLY
    if output_view is OutputView.PRECEDENT_ONLY and source in {SourceGroup.ANY, SourceGroup.MOTIF, SourceGroup.RULES}:
        notes.append("precedent_only output requires literature source execution")

    if run_strategy is RunStrategy.PER_SOURCE:
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
            run_strategy=run_strategy,
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
        run_strategy=run_strategy,
        output_view=output_view,
        notes=tuple(notes),
    )
