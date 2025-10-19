"""
chemtools.formatters - Output formatting for chemistry recommendations.

This package provides modular formatting functions for structuring
recommendation outputs in standardized JSON formats.

Modules:
    base: Core formatting for metadata, input, and detection sections
    normalization: Normalization helpers for chemicals, conditions, recommendations
    rule_output: Rule-based output formatting (SCDB)
    ml_output: ML output formatting and standard builders
    utils: Utility functions for reagent enrichment and condition formatting
"""

# Re-export all public functions for backwards compatibility

# From base.py
from chemtools.formatters.base import (
    format_meta,
    format_input,
    format_detection,
)

# From normalization.py
from chemtools.formatters.normalization import (
    normalize_chemical_entry,
    normalize_condition_value,
    normalize_conditions_block,
    normalize_recommendation_entry,
    normalize_recommendations,
    parse_amount_to_equivalents,
    normalize_rule_string_value,
)

# From rule_output.py
from chemtools.formatters.rule_output import (
    starting_material_entries,
    convert_rule_match_to_recommendations,
)

# From ml_output.py
from chemtools.formatters.ml_output import (
    build_standard_output,
    ensure_standard_output,
    format_ml_output,
    format_rule_output,
    format_fusion_output,
    format_rule_match_result,
)

# From utils.py
from chemtools.formatters.utils import (
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
    
    # Rule output
    "starting_material_entries",
    "convert_rule_match_to_recommendations",
    
    # ML output
    "build_standard_output",
    "ensure_standard_output",
    "format_ml_output",
    "format_rule_output",
    "format_fusion_output",
    "format_rule_match_result",
    
    # Utils
    "enrich_reagent",
    "format_conditions",
    "format_recommendation",
    "parse_condition_string",
]
