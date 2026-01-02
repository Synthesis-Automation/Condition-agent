"""
LangChain integration for ChemTools featurization and analysis.
"""

from .chemtools_wrapper import CHEMTOOLS_TOOLS, get_tool_descriptions, print_tool_summary
from .chemtools_agent import ChemToolsAgent, create_agent, quick_query

__version__ = "2.0.0"

__all__ = [
    "ChemToolsAgent",
    "create_agent",
    "quick_query",
    "CHEMTOOLS_TOOLS",
    "get_tool_descriptions",
    "print_tool_summary",
]
