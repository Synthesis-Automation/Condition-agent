"""Standalone condition-recommendation data and modeling package."""

from .compatibility import CompatibilityAssessment, assess_recipe_compatibility
from .generic_api import (
    GenericConditionRecommender,
    recommend_generic_conditions,
    recommend_indexed_signature,
)
from .label_api import recommend_conditions_from_labels
from .models import (
    AdmissionTier,
    ChemistryStatus,
    ConditionIdentity,
    ConditionStatus,
    GenericConditionRecommendation,
    GenericRecommendationResult,
    IndexEligibility,
    LabelConditionRecommendation,
    LabelRecommendationResult,
    OutcomeStatus,
    ReferenceIdentity,
    RecommendationRecord,
    RecommendationScoreTrace,
    RetrievalLevelTrace,
)
from .rules import RuleRecommendationResult, recommend_rule_conditions
from .rule_review import build_rule_review_rows, export_rule_review_csv

__all__ = [
    "AdmissionTier",
    "ChemistryStatus",
    "ConditionStatus",
    "ConditionIdentity",
    "CompatibilityAssessment",
    "GenericConditionRecommendation",
    "GenericConditionRecommender",
    "GenericRecommendationResult",
    "IndexEligibility",
    "LabelConditionRecommendation",
    "LabelRecommendationResult",
    "OutcomeStatus",
    "ReferenceIdentity",
    "RecommendationRecord",
    "RecommendationScoreTrace",
    "RetrievalLevelTrace",
    "RuleRecommendationResult",
    "build_rule_review_rows",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
    "recommend_rule_conditions",
    "recommend_conditions_from_labels",
    "assess_recipe_compatibility",
    "export_rule_review_csv",
]
