"""
Rule-Based Recommendation System
=================================

Feature-driven condition recommendation for synthetic chemistry reactions.

This module provides a rule-based system that:
1. Analyzes reaction SMILES using the calculable features system
2. Matches features against rule databases (e.g., suzuki.json)
3. Selects appropriate base rules and applies modifiers
4. Generates condition recommendations

Main Components:
- RuleEngine: Core recommendation engine
- RuleDatabase: Rule database loader and validator
- FeatureAnalyzer: Reaction feature extraction
- ConditionRecommendation: Structured recommendation output

Example Usage:
    >>> from chemtools.rule import RuleEngine
    >>> engine = RuleEngine.from_database("suzuki")
    >>> recommendation = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    >>> print(recommendation.to_dict())
"""

from .engine import RuleEngine
from .database import RuleDatabase
from .analyzer import FeatureAnalyzer
from .models import ConditionRecommendation, AppliedRule, AppliedModifier
from .builder import RuleBuilder, ValidationIssue

__all__ = [
    "RuleEngine",
    "RuleDatabase",
    "FeatureAnalyzer",
    "ConditionRecommendation",
    "AppliedRule",
    "AppliedModifier",
    "RuleBuilder",
    "ValidationIssue",
]

__version__ = "1.0.0"
