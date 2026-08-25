"""Deterministic route-refinement contract regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

from core_retrosynthesis import (
    RouteRefinementIntent,
    build_route_refinement_plan,
    collect_route_refinement_issues,
    plan_multistep_routes,
    summarize_route_refinement,
)

from .test_multistep import (
    _LiteratureIndex,
    _candidate,
    _condition_evidence,
    _expander,
)


def _identified_candidate(
    product: str,
    precursors: str,
    identity: str,
    *,
    realization: str,
    score: float,
):
    return replace(
        _candidate(product, precursors, score=score),
        operator_id=f"operator:{identity}",
        disconnection_site_key=f"site:{identity}",
        synthon_signature=f"synthon:{identity}",
        realization_id=realization,
        strategy_id="",
    )


def _condition_evaluator(reaction_smiles: str):
    status = (
        "insufficient_evidence"
        if reaction_smiles.startswith("C.C>>")
        else "recommended_direct"
    )
    return _condition_evidence(reaction_smiles, status)


def test_refinement_excludes_source_strategy_and_finds_improved_route() -> None:
    source_candidate = _identified_candidate(
        "CCCCCCCC",
        "C.C",
        "source",
        realization="realization:source",
        score=0.99,
    )
    alternative = _identified_candidate(
        "CCCCCCCC",
        "N.O",
        "alternative",
        realization="realization:alternative",
        score=0.90,
    )
    expander = _expander({"CCCCCCCC": (source_candidate, alternative)})
    baseline = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        top_k_routes=2,
        condition_evidence_evaluator=_condition_evaluator,
        expander=expander,
    )
    source_route = next(
        route
        for route in baseline.routes
        if route.steps[0].candidate.strategy_id == source_candidate.strategy_id
    )
    issues = collect_route_refinement_issues(source_route)
    condition_issue = next(item for item in issues if item.kind == "condition_gap")
    intent = RouteRefinementIntent(
        source_route_id=source_route.route_id,
        source_step_id=source_route.steps[0].step_id,
        objective="resolve_condition_gap",
        method="alternate_disconnection",
        issue_ids=(condition_issue.issue_id,),
    )
    plan = build_route_refinement_plan(source_route, intent, issues)

    refined = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        top_k_routes=2,
        condition_evidence_evaluator=_condition_evaluator,
        candidate_exclusions=(plan.exclusion,),
        expander=expander,
    )
    outcome = summarize_route_refinement(source_route, intent, refined)

    assert refined.diagnostics.refinement_excluded_candidates == 1
    assert all(
        step.candidate.strategy_id != source_candidate.strategy_id
        for route in refined.routes
        for step in route.steps
    )
    assert outcome.status == "improved_alternative_found"
    assert outcome.source_route_preserved is True
    assert outcome.candidate_assessments[0].relevant_issue_count == 0


def test_alternate_realization_excludes_only_exact_realization() -> None:
    first = _identified_candidate(
        "CCCCCCCC",
        "C.C",
        "shared",
        realization="realization:first",
        score=0.99,
    )
    second = replace(
        first,
        precursor_smiles="N.O",
        proposed_reaction_smiles="N.O>>CCCCCCCC",
        condition_query_reaction_smiles="N.O>>CCCCCCCC",
        realization_id="realization:second",
        score=0.90,
    )
    expander = _expander({"CCCCCCCC": (first, second)})
    baseline = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        top_k_routes=2,
        condition_evidence_evaluator=_condition_evaluator,
        expander=expander,
    )
    source_route = next(
        route
        for route in baseline.routes
        if route.steps[0].candidate.realization_id == "realization:first"
    )
    issue = next(
        item
        for item in collect_route_refinement_issues(source_route)
        if item.kind == "condition_gap"
    )
    intent = RouteRefinementIntent(
        source_route_id=source_route.route_id,
        source_step_id=source_route.steps[0].step_id,
        objective="resolve_condition_gap",
        method="alternate_realization",
        issue_ids=(issue.issue_id,),
    )
    plan = build_route_refinement_plan(
        source_route,
        intent,
        collect_route_refinement_issues(source_route),
    )
    refined = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        candidate_exclusions=(plan.exclusion,),
        expander=expander,
    )

    assert refined.routes[0].steps[0].candidate.strategy_id == first.strategy_id
    assert refined.routes[0].steps[0].candidate.realization_id == "realization:second"


def test_refinement_rejects_objective_issue_mismatch() -> None:
    candidate = _identified_candidate(
        "CCCCCCCC",
        "C.C",
        "source",
        realization="realization:source",
        score=0.99,
    )
    baseline = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=50.0,
        condition_evidence_evaluator=_condition_evaluator,
        expander=_expander({"CCCCCCCC": (candidate,)}),
    )
    route = baseline.routes[0]
    issues = collect_route_refinement_issues(route)
    condition_issue = next(item for item in issues if item.kind == "condition_gap")
    intent = RouteRefinementIntent(
        source_route_id=route.route_id,
        source_step_id=route.steps[0].step_id,
        objective="resolve_compatibility_conflict",
        method="alternate_disconnection",
        issue_ids=(condition_issue.issue_id,),
    )

    with pytest.raises(ValueError, match="not supported by the cited issue"):
        build_route_refinement_plan(route, intent, issues)
