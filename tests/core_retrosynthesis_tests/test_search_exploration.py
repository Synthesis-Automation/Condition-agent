"""Behavioral regressions for widening and dependency-aware guidance."""

from dataclasses import replace
from types import SimpleNamespace as NS

import pytest

from core_retrosynthesis.multistep import (
    _Leaf,
    _RouteState,
    _guided_state_priority,
    StartingMaterialAssessment,
    plan_multistep_routes,
)
from core_retrosynthesis.route_state_learning import LiteratureRouteOrderingGuidance
from core_retrosynthesis.search_exploration import load_search_exploration_policy
from .test_multistep import _candidate, _LiteratureIndex


def _ordering() -> LiteratureRouteOrderingGuidance:
    return LiteratureRouteOrderingGuidance(
        NS(
            train_reverse_operator_transition_support={"A>B": 10, "B>C": 5},
            train_reverse_operator_transition_patent_support={"A>B": 2, "B>C": 2},
        )
    )


def _step(operator: str, product: str, children: tuple[str, ...]) -> NS:
    return NS(
        candidate=NS(operator_id=operator),
        product_node_id=product,
        precursor_node_ids=children,
    )


def test_ordering_uses_dependencies_and_ignores_sibling_traversal() -> None:
    guidance = _ordering()
    a = _step("A", "root", ("left", "right"))
    b = _step("B", "left", ("l1",))
    c = _step("C", "right", ("r1",))
    assert guidance.state_priority(NS(steps=(a, b, c))) == (
        guidance.state_priority(NS(steps=(a, c, b)))
    )
    assert guidance.state_priority(NS(steps=(a, b))) < (0.0,)
    assert guidance.state_priority(NS(steps=(b, c))) == (0.0,)
    # Conflicting producers do not silently establish a dependency.
    duplicate = _step("A", "other", ("left",))
    assert guidance.state_priority(NS(steps=(a, duplicate, b))) == (0.0,)
    assert guidance.state_priority(
        NS(
            steps=(
                _step("A", "", ("",)),
                _step("B", "", ()),
            )
        )
    ) == (0.0,)


def test_guidance_cannot_override_blocked_status_or_cost_band() -> None:
    leaf = StartingMaterialAssessment(
        smiles="CC",
        canonical_smiles="CC",
        depth=1,
        molecular_weight=30,
        terminal=False,
        terminal_reasons=(),
        terminal_evidence="none",
        catalog_role_status="no_match",
    )
    cheap = _RouteState((), (_Leaf(leaf, ()),), 1.0)
    expensive = replace(cheap, cost=100.0)
    blocked = replace(
        cheap, leaves=(_Leaf(replace(leaf, unresolved_reason="no_candidates"), ()),)
    )

    class Guidance:
        def state_priority(self, state):
            return (-1000000 if state.leaves[0].unresolved_reason else 0,)

    def priority(state):
        return _guided_state_priority("CC", state, 3, Guidance())

    assert priority(cheap) < priority(blocked)
    assert priority(cheap) < priority(expensive)


def test_widening_recovers_later_action_without_revalidating_cached_pool() -> None:
    calls = []
    root = "CCCCCCCC"
    candidates = (
        _candidate(root, "CCCCCCC"),
        _candidate(root, "CCCCCC"),
        _candidate(root, "CC.CCO"),
    )

    def expand(product, top_k):
        calls.append((product, top_k))
        return candidates[:top_k] if product == root else ()

    kwargs = dict(
        max_depth=3,
        per_step_top_k=1,
        max_expansions=8,
        molecular_weight_threshold=50,
        expander=expand,
    )
    baseline = plan_multistep_routes(root, object(), _LiteratureIndex(), **kwargs)
    assert not baseline.routes
    calls.clear()
    widened = plan_multistep_routes(
        root,
        object(),
        _LiteratureIndex(),
        widening_factor=3,
        **kwargs,
    )
    assert widened.routes
    assert [x for x in calls if x[0] == root] == [(root, 3)]
    assert widened.diagnostics.widening_revisits == 2
    assert widened.diagnostics.expanded_states <= 8
    assert widened.routes[0].steps[0].candidate.precursor_smiles == "CC.CCO"
    assert dict(widened.routes[0].steps[0].step_cost_components)[
        "candidate_rank_tiebreak"
    ] == pytest.approx(3e-6)
    repeated = plan_multistep_routes(
        root,
        object(),
        _LiteratureIndex(),
        widening_factor=3,
        **kwargs,
    )
    assert widened.to_dict() == repeated.to_dict()


def test_widening_retains_budget_exhaustion_and_rejects_invalid_later_action() -> None:
    root = "CCCCCCCC"
    options = (
        _candidate(root, "CCCCCCC"),
        replace(_candidate(root, "CC.CCO"), forward_validation_status="unresolved"),
    )

    def expand(product, top_k):
        return options[:top_k] if product == root else ()

    kwargs = dict(
        max_depth=3,
        per_step_top_k=1,
        molecular_weight_threshold=50,
        expander=expand,
        widening_factor=2,
    )
    limited = plan_multistep_routes(
        root,
        object(),
        _LiteratureIndex(),
        max_expansions=1,
        **kwargs,
    )
    assert limited.diagnostics.stopped_by_expansion_limit
    assert limited.diagnostics.deferred_widening_states == 1
    assert limited.partial_routes
    assert all(route.steps for route in limited.partial_routes)
    complete = plan_multistep_routes(
        root,
        object(),
        _LiteratureIndex(),
        max_expansions=10,
        **kwargs,
    )
    assert not complete.routes
    assert complete.diagnostics.rejected_invalid_candidates == 1


def test_exploration_policy_and_invalid_widening() -> None:
    assert load_search_exploration_policy().guidance_cost_band == 1.0
    with pytest.raises(ValueError, match="widening factor"):
        plan_multistep_routes("CC", object(), _LiteratureIndex(), widening_factor=0)


def test_plan_routes_cli_passes_widening_to_canonical_planner(
    monkeypatch, capsys
) -> None:
    from contextlib import nullcontext
    from core_retrosynthesis import cli

    received = {}

    def planner(*args, **kwargs):
        received.update(kwargs)
        return NS(to_dict=lambda: {"ok": True})

    monkeypatch.setattr(cli, "load_generic_library", lambda path: object())
    monkeypatch.setattr(cli, "open_stock_lookup", lambda path: nullcontext(object()))
    monkeypatch.setattr(cli, "plan_multistep_routes", planner)
    assert (
        cli.main(
            [
                "plan-routes",
                "library",
                "stock",
                "CC",
                "--max-depth",
                "6",
                "--widening-factor",
                "3",
            ]
        )
        == 0
    )
    assert received["widening_factor"] == 3
    assert received["max_depth"] == 6
