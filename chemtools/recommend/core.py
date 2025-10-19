"""
Core recommendation engine using DRFP-based reaction similarity.

REFACTORED: The original 1,334-line monolithic file has been split into
logical modules in chemtools/recommend/modules/:
    - recommender.py: Main DRFP-based recommendation engine
    - structured.py: API-friendly wrapper with enhanced metadata
    - fusion_adapter.py: Convert fusion recommender output to core format
    - precedent_builder.py: Build precedent summaries and statistics
    - output_builder.py: Format multi-variant outputs

This module maintains backwards compatibility by re-exporting all public
functions from the modules package. For new code, prefer importing from
modules directly.

Key Philosophy:
- Use reaction-level similarity (DRFP) as the primary method
- Avoid complex substrate analysis and family-specific featurization
- Keep the code general and maintainable
- Delegate to precedent.knn() for similarity search

Public Functions:
    - recommend_from_reaction(): Main recommendation function
    - recommend_conditions_structured(): Structured output with multiple variants
"""

from __future__ import annotations

# Re-export all public functions from modules package
from .modules import (
    recommend_from_reaction,
    recommend_conditions_structured,
)

# Also re-export internal functions for any code that might use them
# (though these should be considered private)
from .modules.fusion_adapter import (
    convert_fusion_to_core_format as _convert_fusion_to_core_format,
    build_formatted_output_from_fusion as _build_formatted_output_from_fusion,
)
from .modules.precedent_builder import (
    build_precedent_details as _build_precedent_details,
    calculate_average_yield as _calculate_average_yield,
    calculate_yield_range as _calculate_yield_range,
    calculate_temp_range as _calculate_temp_range,
    calculate_time_range as _calculate_time_range,
)
from .modules.output_builder import (
    build_formatted_output as _build_formatted_output,
)

__all__ = [
    # Main public API
    "recommend_from_reaction",
    "recommend_conditions_structured",
    
    # Internal functions (for backwards compatibility only)
    "_convert_fusion_to_core_format",
    "_build_formatted_output_from_fusion",
    "_build_precedent_details",
    "_calculate_average_yield",
    "_calculate_yield_range",
    "_calculate_temp_range",
    "_calculate_time_range",
    "_build_formatted_output",
]

# Note: The actual implementations have been moved to chemtools/recommend/modules/
# This file now serves as a compatibility layer.
#
# For new code, prefer:
#   from chemtools.recommend.modules.recommender import recommend_from_reaction
#   from chemtools.recommend.modules.structured import recommend_conditions_structured
#
# Or simply:
#   from chemtools.recommend import recommend_from_reaction
