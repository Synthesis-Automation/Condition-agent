"""
Reaction condition recommendation package.

The unified recommender combines reaction datasets, literature protocols,
and HTE screens into a single DRFP + feature-tag similarity engine.
"""

from .unified import (
    UnifiedRecommender,
    RecommendationResult,
    recommend_from_reaction,
    recommend_conditions_structured,
)
from .index_builder import build_unified_recommendation_index

__all__ = [
    "UnifiedRecommender",
    "RecommendationResult",
    "recommend_from_reaction",
    "recommend_conditions_structured",
    "build_unified_recommendation_index",
]
