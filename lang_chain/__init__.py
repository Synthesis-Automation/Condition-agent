"""
LangChain integration for ChemTools.

This package provides LangChain/LangGraph wrappers for ChemTools functions,
enabling AI-powered chemistry analysis through a ReAct agent.

Quick Start:
    # Interactive CLI
    python -m lang_chain.chemtools_cli
    
    # Python API
    from lang_chain.chemtools_agent import quick_query
    result = quick_query("Recommend conditions for Suzuki coupling")

Main Modules:
    - chemtools_wrapper: LangChain @tool decorators for chemtools functions
    - chemtools_agent: ReAct agent with LangGraph
    - chemtools_cli: Interactive command-line interface

Available Tools:
    - normalize_smiles_tool
    - normalize_reaction_tool
    - detect_reaction_family_tool
    - classify_reactant_tool
    - get_functional_groups_tool
    - recommend_conditions_tool
    - search_precedents_tool
    - find_reagent_tool

See README.md for complete documentation.
"""

from .chemtools_wrapper import (
    CHEMTOOLS_TOOLS,
    get_tool_descriptions,
    print_tool_summary,
    clear_recommendation_cache,
    recommendation_cache_stats,
)
from .chemtools_agent import ChemToolsAgent, create_agent, quick_query

__version__ = "1.0.0"

__all__ = [
    # Agent
    "ChemToolsAgent",
    "create_agent",
    "quick_query",
    
    # Tools
    "CHEMTOOLS_TOOLS",
    "get_tool_descriptions",
    "print_tool_summary",
    "clear_recommendation_cache",
    "recommendation_cache_stats",
]
