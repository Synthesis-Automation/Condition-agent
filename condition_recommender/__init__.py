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
    ChemistRankingPreferences,
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
from .discovery import ReactionDiscoveryExplorer, load_discovery_rules
from .discovery_models import (
    DiscoveryScoreTrace,
    ReactionDiscoveryHit,
    ReactionDiscoveryResult,
)
from .ranking_preferences import (
    available_ranking_profiles,
    load_chemist_ranking_profiles,
    resolve_ranking_preferences,
)

__all__ = [
    "AdmissionTier",
    "ChemistRankingPreferences",
    "ChemistryStatus",
    "ConditionStatus",
    "ConditionStageStatus",
    "ConditionIdentity",
    "FragmentSourceSupport",
    "CompatibilityAssessment",
    "CoreEligibilityAssessment",
    "CoreEligibility",
    "DiscoveryScoreTrace",
    "GenericConditionRecommendation",
    "GenericConditionRecommender",
    "GenericRecommendationResult",
    "ReactionDiscoveryExplorer",
    "ReactionDiscoveryHit",
    "ReactionDiscoveryResult",
    "IndexEligibility",
    "OutcomeStatus",
    "PrecedentIndexScope",
    "PrecedentTier",
    "ReferenceIdentity",
    "RecommendationRecord",
    "RecommendationScoreTrace",
    "RetrievalLevelTrace",
    "available_ranking_profiles",
    "evaluate_reaction_core_index",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
    "load_chemist_ranking_profiles",
    "load_discovery_rules",
    "resolve_ranking_preferences",
    "assess_recipe_compatibility",
    "assess_core_eligibility",
]
