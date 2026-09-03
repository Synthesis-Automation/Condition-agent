"""Declarative policy for fair continuation of strategic route entries."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Mapping


PORTFOLIO_CONTINUATION_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "portfolio_continuation.v1.json"
)
PORTFOLIO_CONTINUATION_VERSION = "portfolio_continuation.v1"


@dataclass(frozen=True)
class PortfolioContinuationPolicy:
    """Validated bounded-search quotas for first-action route lanes."""

    definition_id: str
    schema_version: str
    minimum_expansions_per_first_action: int
    maximum_active_first_actions: int
    scheduling_mode: str
    ranking_influence: str

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible policy view."""

        return {
            "definition_id": self.definition_id,
            "schema_version": self.schema_version,
            "minimum_expansions_per_first_action": (
                self.minimum_expansions_per_first_action
            ),
            "maximum_active_first_actions": self.maximum_active_first_actions,
            "scheduling_mode": self.scheduling_mode,
            "ranking_influence": self.ranking_influence,
        }


def validate_portfolio_continuation_policy(value: Mapping[str, Any]) -> None:
    """Reject unsafe, unbounded, or ranking-active continuation policies."""

    if value.get("definition_id") != PORTFOLIO_CONTINUATION_VERSION:
        raise ValueError("unexpected portfolio-continuation definition ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported portfolio-continuation schema")
    minimum = value.get("minimum_expansions_per_first_action")
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
        raise ValueError("minimum continuation expansions must be positive")
    maximum = value.get("maximum_active_first_actions")
    if isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1:
        raise ValueError("maximum active first actions must be positive")
    if maximum > 16:
        raise ValueError("maximum active first actions exceeds review bound")
    if value.get("scheduling_mode") != (
        "least_expanded_lane_then_state_priority"
    ):
        raise ValueError("unsupported portfolio-continuation scheduling mode")
    if value.get("ranking_influence") != "review_only_opt_in":
        raise ValueError("portfolio continuation must remain review-only")


@lru_cache(maxsize=1)
def load_portfolio_continuation_policy() -> PortfolioContinuationPolicy:
    """Load the canonical first-action continuation policy."""

    value = json.loads(
        PORTFOLIO_CONTINUATION_DEFINITION_PATH.read_text(encoding="utf-8")
    )
    validate_portfolio_continuation_policy(value)
    return PortfolioContinuationPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        minimum_expansions_per_first_action=int(
            value["minimum_expansions_per_first_action"]
        ),
        maximum_active_first_actions=int(value["maximum_active_first_actions"]),
        scheduling_mode=str(value["scheduling_mode"]),
        ranking_influence=str(value["ranking_influence"]),
    )


__all__ = [
    "PORTFOLIO_CONTINUATION_DEFINITION_PATH",
    "PORTFOLIO_CONTINUATION_VERSION",
    "PortfolioContinuationPolicy",
    "load_portfolio_continuation_policy",
    "validate_portfolio_continuation_policy",
]
