"""Standalone condition-recommendation data and modeling package."""

from .compatibility import CompatibilityAssessment, assess_recipe_compatibility
from .core_evaluation import evaluate_reaction_core_index
from .generic_api import (
    GenericConditionRecommender,
    recommend_generic_conditions,
    recommend_indexed_signature,
)
from .models import (
    AdmissionTier,
    ChemistryStatus,
    ConditionIdentity,
    ConditionStageStatus,
    ConditionStatus,
    FragmentSourceSupport,
    GenericConditionRecommendation,
    GenericRecommendationResult,
    IndexEligibility,
    OutcomeStatus,
    ReferenceIdentity,
    RecommendationRecord,
    RecommendationScoreTrace,
    RetrievalLevelTrace,
)

__all__ = [
    "AdmissionTier",
    "ChemistryStatus",
    "ConditionStatus",
    "ConditionStageStatus",
    "ConditionIdentity",
    "FragmentSourceSupport",
    "CompatibilityAssessment",
    "GenericConditionRecommendation",
    "GenericConditionRecommender",
    "GenericRecommendationResult",
    "IndexEligibility",
    "OutcomeStatus",
    "ReferenceIdentity",
    "RecommendationRecord",
    "RecommendationScoreTrace",
    "RetrievalLevelTrace",
    "evaluate_reaction_core_index",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
    "assess_recipe_compatibility",
]
