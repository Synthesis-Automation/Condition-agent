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
    ReactionCompletionOption,
    ReactionCompletionProposal,
    ReactionCompletionRequirement,
    ReactionCompletionSelection,
    ReferenceIdentity,
    RecommendationRecord,
    RecommendationScoreTrace,
    RetrievalLevelTrace,
)
from .ranking_preferences import (
    available_ranking_profiles,
    load_chemist_ranking_profiles,
    resolve_ranking_preferences,
)
from .reaction_completion import (
    build_completion_selection,
    build_completed_reaction_smiles,
    propose_reaction_completion,
)
from .recipe_assessment import assess_reaction_recipe

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
    "GenericConditionRecommendation",
    "GenericConditionRecommender",
    "GenericRecommendationResult",
    "IndexEligibility",
    "OutcomeStatus",
    "PrecedentIndexScope",
    "PrecedentTier",
    "ReactionCompletionOption",
    "ReactionCompletionProposal",
    "ReactionCompletionRequirement",
    "ReactionCompletionSelection",
    "assess_reaction_recipe",
    "ReferenceIdentity",
    "RecommendationRecord",
    "RecommendationScoreTrace",
    "RetrievalLevelTrace",
    "available_ranking_profiles",
    "evaluate_reaction_core_index",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
    "load_chemist_ranking_profiles",
    "resolve_ranking_preferences",
    "build_completion_selection",
    "build_completed_reaction_smiles",
    "propose_reaction_completion",
    "assess_recipe_compatibility",
    "assess_core_eligibility",
]
