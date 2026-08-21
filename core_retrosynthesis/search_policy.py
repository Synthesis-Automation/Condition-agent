"""Typed bounds for one authorized multistep search expansion."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict


ROUTE_SEARCH_POLICY_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class RouteSearchPolicy:
    """Initial search settings and user-authorized maximum settings."""

    initial_max_depth: int = 2
    initial_beam_width: int = 10
    initial_max_expansions: int = 6
    maximum_max_depth: int = 3
    maximum_beam_width: int = 20
    maximum_max_expansions: int = 12
    schema_version: str = ROUTE_SEARCH_POLICY_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.initial_max_depth not in {2, 3} or self.maximum_max_depth not in {2, 3}:
            raise ValueError("route depth must be two or three")
        if self.initial_max_depth > self.maximum_max_depth:
            raise ValueError("initial route depth exceeds authorized maximum")
        for initial, maximum, name in (
            (self.initial_beam_width, self.maximum_beam_width, "beam width"),
            (self.initial_max_expansions, self.maximum_max_expansions, "expansions"),
        ):
            if initial < 1 or maximum < initial:
                raise ValueError(f"invalid initial/maximum route {name}")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RouteSearchPolicyDelta:
    """A monotonic bounded enlargement; never a route edit."""

    max_depth_delta: int = 0
    beam_width_delta: int = 0
    max_expansions_delta: int = 0
    schema_version: str = ROUTE_SEARCH_POLICY_SCHEMA_VERSION

    def __post_init__(self) -> None:
        values = (
            self.max_depth_delta,
            self.beam_width_delta,
            self.max_expansions_delta,
        )
        if any(value < 0 for value in values):
            raise ValueError("route search expansion deltas cannot be negative")
        if not any(values):
            raise ValueError("route search expansion delta must add capacity")


def apply_route_search_delta(
    policy: RouteSearchPolicy,
    delta: RouteSearchPolicyDelta,
    *,
    current_max_depth: int,
    current_beam_width: int,
    current_max_expansions: int,
) -> tuple[int, int, int]:
    """Return enlarged settings or reject authority/bound violations."""

    values = (
        current_max_depth + delta.max_depth_delta,
        current_beam_width + delta.beam_width_delta,
        current_max_expansions + delta.max_expansions_delta,
    )
    maxima = (
        policy.maximum_max_depth,
        policy.maximum_beam_width,
        policy.maximum_max_expansions,
    )
    if any(value > maximum for value, maximum in zip(values, maxima)):
        raise ValueError("route search expansion exceeds user-authorized maximum")
    if values[0] not in {2, 3}:
        raise ValueError("expanded route depth must be two or three")
    return values
