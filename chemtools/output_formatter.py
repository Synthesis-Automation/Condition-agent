"""
Enhanced output formatter for reaction condition recommendations.

This module maintains backwards compatibility by re-exporting all formatting
functions from the chemtools.formatters package.

REFACTORED: The original 1,398-line monolithic file has been split into
logical modules in chemtools/formatters/:
    - base.py: Core formatting (metadata, input, detection)
    - normalization.py: Normalization helpers for chemicals and conditions
    - ml_output.py: ML output formatting and standard builders
    - utils.py: Reagent enrichment and condition formatting utilities

For new code, prefer importing directly from chemtools.formatters submodules.
This file provides a compatibility layer for existing code.
"""

# Re-export all public functions from formatters package
from chemtools.formatters import (
    # Base formatting (metadata, input, detection)
    format_meta,
    format_input,
    format_detection,
    
    # Normalization helpers
    normalize_chemical_entry,
    normalize_condition_value,
    normalize_conditions_block,
    normalize_recommendation_entry,
    normalize_recommendations,
    parse_amount_to_equivalents,
    normalize_rule_string_value,
    
    # ML output and standard builders
    build_standard_output,
    ensure_standard_output,
    format_ml_output,
    format_fusion_output,
    
    # Utility functions
    enrich_reagent,
    format_conditions,
    format_recommendation,
    parse_condition_string,
)

__all__ = [
    # Base
    "format_meta",
    "format_input",
    "format_detection",
    
    # Normalization
    "normalize_chemical_entry",
    "normalize_condition_value",
    "normalize_conditions_block",
    "normalize_recommendation_entry",
    "normalize_recommendations",
    "parse_amount_to_equivalents",
    "normalize_rule_string_value",
    
    # ML output
    "build_standard_output",
    "ensure_standard_output",
    "format_ml_output",
    "format_fusion_output",
    
    # Utils
    "enrich_reagent",
    "format_conditions",
    "format_recommendation",
    "parse_condition_string",
]
