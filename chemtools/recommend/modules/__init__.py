"""
chemtools.recommend.modules - Modular recommendation components.

This package contains the refactored recommendation engine split into
logical modules for better maintainability and testability.

Modules:
    recommender: Main DRFP-based recommendation engine
    structured: API-friendly wrapper with enhanced metadata
    fusion_adapter: Convert fusion recommender output to core format
    precedent_builder: Build precedent summaries and statistics
    output_builder: Format multi-variant outputs
"""

from .recommender import recommend_from_reaction
from .structured import recommend_conditions_structured
from .fusion_adapter import (
    convert_fusion_to_core_format,
    build_formatted_output_from_fusion,
)
from .precedent_builder import (
    build_precedent_details,
    calculate_average_yield,
    calculate_yield_range,
    calculate_temp_range,
    calculate_time_range,
)
from .output_builder import build_formatted_output

__all__ = [
    # Main public functions
    "recommend_from_reaction",
    "recommend_conditions_structured",
    
    # Fusion adapter
    "convert_fusion_to_core_format",
    "build_formatted_output_from_fusion",
    
    # Precedent builders
    "build_precedent_details",
    "calculate_average_yield",
    "calculate_yield_range",
    "calculate_temp_range",
    "calculate_time_range",
    
    # Output builder
    "build_formatted_output",
]
