"""
Reaction condition recommendation package.

This package provides DRFP-based reaction similarity search for recommending
optimal reaction conditions. The core philosophy is to use reaction-level
fingerprinting (DRFP) rather than complex substrate featurization for better
generality and maintainability.

Public API:
    - recommend_from_reaction(): Main recommendation function
    - recommend_conditions_structured(): Structured output for API
    - design_plate_from_reaction(): Experimental plate design (from plate_design module)

Example:
    >>> from chemtools.recommend import recommend_from_reaction
    >>> results = recommend_from_reaction(
    ...     reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    ...     k=50
    ... )
    >>> print(results['recommendation'])
    {'core': 'Pd/XPhos', 'base_uid': '7758-29-4', 'solvent_uid': '1120-21-4', ...}
"""

from .core import recommend_from_reaction, recommend_conditions_structured
from .hte_simple import HTESimpleConditionRecommender, RecommendationOptions
from .plate_design import design_plate_from_reaction

__all__ = [
    "recommend_from_reaction",
    "recommend_conditions_structured",
    "HTESimpleConditionRecommender",
    "RecommendationOptions",
    "design_plate_from_reaction",
]
