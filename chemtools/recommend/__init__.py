"""
Reaction condition recommendation package.

Primary backend: chemtools.recommend (High-Throughput Experimentation).
"""

from __future__ import annotations

from typing import Any

from .hte_adapter import recommend_from_reaction, recommend_conditions_structured

_RECOMMENDER_EXPORTS = {
    "HTERecommender",
    "HTERecommendationResult",
    "ConditionRecommendation",
    "format_recommendation",
    "format_result",
}

_ANALYTICS_EXPORTS = {"HTEAnalytics"}


def __getattr__(name: str) -> Any:
    if name in _RECOMMENDER_EXPORTS:
        from . import recommender as _recommender
        return getattr(_recommender, name)
    if name in _ANALYTICS_EXPORTS:
        from . import analytics as _analytics
        return getattr(_analytics, name)
    raise AttributeError(f"module 'chemtools.recommend' has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(
        list(globals().keys())
        + list(_RECOMMENDER_EXPORTS)
        + list(_ANALYTICS_EXPORTS)
    )


__all__ = [
    "recommend_from_reaction",
    "recommend_conditions_structured",
    "HTERecommender",
    "HTERecommendationResult",
    "ConditionRecommendation",
    "format_recommendation",
    "format_result",
    "HTEAnalytics",
]
