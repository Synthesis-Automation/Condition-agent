"""Standalone condition-recommendation data and modeling package."""

from .compatibility import CompatibilityAssessment, assess_recipe_compatibility
from .core_evaluation import evaluate_reaction_core_index
from .core_eligibility import (
    CoreEligibilityAssessment,
    assess_core_eligibility,
)
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
    CoreEligibility,
    FragmentSourceSupport,
    GenericConditionRecommendation,
    GenericRecommendationResult,
    IndexEligibility,
    OutcomeStatus,
    PrecedentIndexScope,
    PrecedentTier,
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
    "CoreEligibilityAssessment",
    "CoreEligibility",
    "GenericConditionRecommendation",
    "GenericConditionRecommender",
    "GenericRecommendationResult",
    "IndexEligibility",
    "OutcomeStatus",
    "PrecedentIndexScope",
    "PrecedentTier",
    "ReferenceIdentity",
    "RecommendationRecord",
    "RecommendationScoreTrace",
    "RetrievalLevelTrace",
    "evaluate_reaction_core_index",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
    "assess_recipe_compatibility",
    "assess_core_eligibility",
]
