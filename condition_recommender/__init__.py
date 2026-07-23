"""Standalone condition-recommendation data and modeling package."""

from .api import recommend_conditions
from .compatibility import CompatibilityAssessment, assess_recipe_compatibility
from .generic_api import recommend_generic_conditions, recommend_indexed_signature
from .label_api import recommend_conditions_from_labels
from .models import (
    AdmissionTier,
    ConditionIdentity,
    ConditionRecommendation,
    GenericConditionRecommendation,
    GenericRecommendationResult,
    LabelConditionRecommendation,
    LabelRecommendationResult,
    RecommendationRecord,
    RecommendationResult,
)
from .rules import RuleRecommendationResult, recommend_rule_conditions
from .rule_review import build_rule_review_rows, export_rule_review_csv

__all__ = [
    "AdmissionTier",
    "ConditionIdentity",
    "ConditionRecommendation",
    "CompatibilityAssessment",
    "GenericConditionRecommendation",
    "GenericRecommendationResult",
    "LabelConditionRecommendation",
    "LabelRecommendationResult",
    "RecommendationRecord",
    "RecommendationResult",
    "RuleRecommendationResult",
    "build_rule_review_rows",
    "recommend_conditions",
    "recommend_generic_conditions",
    "recommend_indexed_signature",
    "recommend_rule_conditions",
    "recommend_conditions_from_labels",
    "assess_recipe_compatibility",
    "export_rule_review_csv",
]
