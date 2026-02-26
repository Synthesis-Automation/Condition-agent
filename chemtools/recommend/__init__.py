"""
Reaction condition recommendation package.

Primary backend: chemtools.recommend (High-Throughput Experimentation).
"""

from __future__ import annotations

from typing import Any

from .hte_adapter import recommend_from_reaction, recommend_conditions_structured
from .api import recommend, warm_recommendation_cache
from .data_manager import RecommendationDataManager
from .models import (
    OutputView,
    QueryAnalysis,
    RecommendationStrategy,
    RecommendationRequest,
    RecommendationRunResult,
    RunStrategy,
    SourceGroup,
    SourcePlan,
)
from .planner import plan_sources
from .query_analysis import analyze_recommendation_query

_RECOMMENDER_EXPORTS = {
    "HTERecommender",
    "HTERecommendationResult",
    "ConditionRecommendation",
    "warm_hte_cache",
    "format_recommendation",
    "format_result",
}

_ANALYTICS_EXPORTS = {"HTEAnalytics"}
_API_EXPORTS = {
    "recommend",
    "warm_recommendation_cache",
    "RecommendationDataManager",
    "RecommendationRequest",
    "QueryAnalysis",
    "SourcePlan",
    "RecommendationRunResult",
    "RecommendationStrategy",
    "SourceGroup",
    "RunStrategy",
    "OutputView",
    "plan_sources",
    "analyze_recommendation_query",
}


def __getattr__(name: str) -> Any:
    if name in _API_EXPORTS:
        return globals()[name]
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
        + list(_API_EXPORTS)
        + list(_RECOMMENDER_EXPORTS)
        + list(_ANALYTICS_EXPORTS)
    )


__all__ = [
    "recommend_from_reaction",
    "recommend_conditions_structured",
    "HTERecommender",
    "HTERecommendationResult",
    "ConditionRecommendation",
    "warm_hte_cache",
    "format_recommendation",
    "format_result",
    "HTEAnalytics",
    "recommend",
    "warm_recommendation_cache",
    "RecommendationDataManager",
    "RecommendationRequest",
    "QueryAnalysis",
    "SourcePlan",
    "RecommendationRunResult",
    "RecommendationStrategy",
    "SourceGroup",
    "RunStrategy",
    "OutputView",
    "plan_sources",
    "analyze_recommendation_query",
]
