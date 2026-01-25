"""
Reaction condition recommendation package.

Primary backend: chemtools.HTE (High-Throughput Experimentation).
"""

from __future__ import annotations

from typing import Any

from .hte_adapter import recommend_from_reaction, recommend_conditions_structured

_HTE_EXPORTS = {
    "HTERecommender",
    "HTERecommendationResult",
    "ConditionRecommendation",
    "format_recommendation",
    "format_result",
    "HTEAnalytics",
}


def __getattr__(name: str) -> Any:
    if name in _HTE_EXPORTS:
        from chemtools import HTE as _hte
        return getattr(_hte, name)
    raise AttributeError(f"module 'chemtools.recommend' has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(list(globals().keys()) + list(_HTE_EXPORTS))


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
