"""Review-only latent-state action portfolio regressions."""

from __future__ import annotations

from dataclasses import replace
import json

import pytest

from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    LATENT_STATE_ROUTE_SEARCH_DEFINITION_PATH,
    LatentStateActionSelector,
    PORTFOLIO_CONTINUATION_DEFINITION_PATH,
    PartitionRealizationConfig,
    classify_route_action,
    load_latent_state_route_search_policy,
    load_portfolio_continuation_policy,
    plan_multistep_routes,
    select_latent_state_action_portfolio,
    validate_latent_state_route_search_policy,
    validate_portfolio_continuation_policy,
)

from .test_multistep import _LiteratureIndex, _candidate


def _mapped_candidate(
    reaction_smiles: str,
    *,
    strategic_class: str = "fragmentation",
) -> GenericDisconnectionCandidate:
    precursors, product = reaction_smiles.split(">>", 1)
    return replace(
        _candidate(product, precursors),
        proposed_reaction_smiles=reaction_smiles,
        condition_query_reaction_smiles=reaction_smiles,
        strategic_class=strategic_class,
    )


CONVERGENT = _mapped_candidate(
    "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]"
)
SINGLE_CARRIER = _mapped_candidate(
    "[CH3:1][CH2:2][CH2:3][OH:4].[CH3:5][Br:6]>>"
    "[CH3:1][CH2:2][CH2:3][O:4][CH3:5]"
)
UNARY = _mapped_candidate(
    "[CH3:1][CH:2]=[O:3]>>[CH3:1][CH2:2][OH:3]"
)
RING = _mapped_candidate(
    "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]",
    strategic_class="ring_construction",
)


def test_latent_state_policy_is_validated_and_review_only() -> None:
    raw = json.loads(
        LATENT_STATE_ROUTE_SEARCH_DEFINITION_PATH.read_text(encoding="utf-8")
    )

    validate_latent_state_route_search_policy(raw)
    assert load_latent_state_route_search_policy().ranking_influence == (
        "review_only_opt_in"
    )
    raw["route_portfolio"]["ranking_influence"] = "default_ranking"
    try:
        validate_latent_state_route_search_policy(raw)
    except ValueError as error:
        assert "review-only" in str(error)
    else:
        raise AssertionError("ranking-active policy must be rejected")


def test_portfolio_continuation_policy_is_bounded_and_review_only() -> None:
    raw = json.loads(
        PORTFOLIO_CONTINUATION_DEFINITION_PATH.read_text(encoding="utf-8")
    )

    validate_portfolio_continuation_policy(raw)
    policy = load_portfolio_continuation_policy()
    assert policy.minimum_expansions_per_first_action == 2
    assert policy.maximum_active_first_actions == 4
    raw["maximum_active_first_actions"] = 17
    with pytest.raises(ValueError, match="review bound"):
        validate_portfolio_continuation_policy(raw)


def test_route_actions_are_classified_from_mapped_target_atom_distribution() -> None:
    assert classify_route_action(CONVERGENT).action_class == "convergent_assembly"
    assert classify_route_action(SINGLE_CARRIER).action_class == (
        "single_carrier_installation"
    )
    assert classify_route_action(UNARY).action_class == "unary_state_change"
    assert classify_route_action(RING).action_class == "ring_reorganization"
    evidence = classify_route_action(SINGLE_CARRIER)
    assert evidence.largest_target_atom_fraction == 0.8
    assert evidence.target_atom_coverage == 1.0


def test_candidate_portfolio_preserves_strategically_distinct_validated_actions() -> None:
    candidates = (
        CONVERGENT,
        SINGLE_CARRIER,
        replace(SINGLE_CARRIER, template_id="single-carrier:2"),
        replace(CONVERGENT, template_id="convergent:2"),
        RING,
        UNARY,
    )

    selected = select_latent_state_action_portfolio(candidates, top_k=4)

    assert selected[0] is CONVERGENT
    selected_classes = tuple(
        classify_route_action(item).action_class for item in selected
    )
    assert set(selected_classes) == {
        "convergent_assembly",
        "single_carrier_installation",
        "ring_reorganization",
    }
    assert selected_classes.count("single_carrier_installation") == 2


def test_multistep_selector_requests_wider_pool_then_keeps_bounded_actions() -> None:
    requested = []
    first = _candidate("CCCCCCCC", "CC.CC")
    second = _candidate("CCCCCCCC", "C.C")

    class _Selector:
        definition_id = "test_selector.v1"
        candidate_pool_multiplier = 3
        continuation_definition_id = "test_continuation.v1"
        minimum_expansions_per_first_action = 1
        maximum_active_first_actions = 2

        def __call__(self, product_smiles, candidates, top_k):
            assert product_smiles == "CCCCCCCC"
            assert top_k == 5
            return candidates[-1:]

        def state_diversity_key(self, state):
            return "selected"

        def continuation_lane_key(self, state):
            return state.steps[0].step_id if state.steps else "root"

        def continuation_lane_class(self, state):
            return "selected"

        def select_routes(self, routes, limit):
            return routes[:limit]

    def expand(product: str, top_k: int):
        requested.append((product, top_k))
        return (first, second)[:top_k]

    result = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        molecular_weight_threshold=50.0,
        per_step_top_k=5,
        max_candidates_to_validate=20,
        expander=expand,
        route_action_selector=_Selector(),
    )

    assert requested == [("CCCCCCCC", 15)]
    assert result.routes[0].steps[0].candidate is second
    assert result.diagnostics.route_action_pool_candidates == 2
    assert result.diagnostics.route_action_selected_candidates == 1
    assert result.diagnostics.route_action_selector_definition_id == (
        "test_selector.v1"
    )


def test_partition_config_keeps_latent_portfolio_opt_in_and_backward_readable() -> None:
    config = PartitionRealizationConfig()

    assert config.use_latent_state_portfolio is False
    serialized = config.to_dict()
    serialized.pop("use_latent_state_portfolio")
    assert PartitionRealizationConfig.from_dict(serialized) == config
    assert LatentStateActionSelector().candidate_pool_multiplier == 4


def test_multistep_selector_cannot_inject_an_unvalidated_action() -> None:
    candidate = _candidate("CCCCCCCC", "C.C")

    class _InjectingSelector:
        definition_id = "invalid_selector.v1"
        candidate_pool_multiplier = 1
        continuation_definition_id = "test_continuation.v1"
        minimum_expansions_per_first_action = 1
        maximum_active_first_actions = 1

        def __call__(self, product_smiles, candidates, top_k):
            del product_smiles, top_k
            return (replace(candidates[0], template_id="injected"),)

        def state_diversity_key(self, state):
            return "injected"

        def continuation_lane_key(self, state):
            return state.steps[0].step_id if state.steps else "root"

        def continuation_lane_class(self, state):
            return "injected"

        def select_routes(self, routes, limit):
            return routes[:limit]

    with pytest.raises(ValueError, match="validated candidate pool"):
        plan_multistep_routes(
            "CCCCCCCC",
            object(),
            _LiteratureIndex(),
            molecular_weight_threshold=50.0,
            expander=lambda product, top_k: (candidate,),
            route_action_selector=_InjectingSelector(),
        )


def test_portfolio_continuation_expands_each_retained_first_action_lane() -> None:
    first_a = _candidate("CCCCCCCCCCCC", "C.CCCCCCCCCC")
    first_b = _candidate("CCCCCCCCCCCC", "C.CCCCCCCCN")
    expansions = {
        "CCCCCCCCCCCC": (first_a, first_b),
        "CCCCCCCCCC": (_candidate("CCCCCCCCCC", "C.CCCCCCCC"),),
        "CCCCCCCCN": (_candidate("CCCCCCCCN", "C.CCCCCCN"),),
    }
    calls = []

    def expand(product: str, top_k: int):
        calls.append(product)
        return expansions.get(product, ())[:top_k]

    result = plan_multistep_routes(
        "CCCCCCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=3,
        molecular_weight_threshold=50.0,
        per_step_top_k=2,
        max_expansions=5,
        expander=expand,
        route_action_selector=LatentStateActionSelector(),
    )

    assert calls[0] == "CCCCCCCCCCCC"
    assert calls.count("CCCCCCCCCC") == 1
    assert calls.count("CCCCCCCCN") == 1
    assert result.diagnostics.continuation_active_lane_count == 2
    assert result.diagnostics.continuation_quota_selected_states == 4
    assert result.diagnostics.continuation_lanes_reaching_minimum == 2
    assert {
        count
        for _, _, count in result.diagnostics.continuation_lane_expansions
    } == {2}
