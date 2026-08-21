"""Bounded route-search expansion policy tests."""

from __future__ import annotations

import pytest

from core_retrosynthesis import (
    RouteSearchPolicy,
    RouteSearchPolicyDelta,
    apply_route_search_delta,
)


def test_route_search_delta_is_monotonic_and_within_authorized_maxima() -> None:
    values = apply_route_search_delta(
        RouteSearchPolicy(),
        RouteSearchPolicyDelta(
            max_depth_delta=1,
            beam_width_delta=5,
            max_expansions_delta=4,
        ),
        current_max_depth=2,
        current_beam_width=10,
        current_max_expansions=6,
    )

    assert values == (3, 15, 10)


def test_route_search_delta_cannot_exceed_original_authority() -> None:
    with pytest.raises(ValueError, match="authorized maximum"):
        apply_route_search_delta(
            RouteSearchPolicy(),
            RouteSearchPolicyDelta(beam_width_delta=11),
            current_max_depth=2,
            current_beam_width=10,
            current_max_expansions=6,
        )


def test_route_search_delta_cannot_reduce_or_repeat_same_search() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        RouteSearchPolicyDelta(max_expansions_delta=-1)
    with pytest.raises(ValueError, match="must add capacity"):
        RouteSearchPolicyDelta()
