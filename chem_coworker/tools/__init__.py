"""
ToolRegistry — the self-describing tool system for ChemCoworker.

How it works:
  - Each tool module (chemistry.py, taxonomy.py, etc.) defines ToolPlugin objects
    and calls REGISTRY.register() at import time
  - The registry builds a dependency graph and can sort tools into parallel
    execution groups via topological sort
  - auto-generates the tool description block for REASON_PROMPT
  - The ToolExecutor uses get_callables() to run tools directly (no LangChain overhead)

Adding a new tool:
  1. Create chem_coworker/tools/my_new_tools.py
  2. Define ToolPlugin objects with fn=, prerequisites=, etc.
  3. Call REGISTRY.register(plugin) for each
  4. Import the module below so it runs at startup
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional

from ._base import ToolPlugin


class ToolRegistry:
    """
    Self-describing registry of ChemCoworker tools.

    Supports:
    - Topological sort of tool names into parallel execution groups
    - Auto-generation of tool descriptions for system prompts
    - Direct callable access (no LangChain overhead for tool execution)
    """

    def __init__(self) -> None:
        self._plugins: Dict[str, ToolPlugin] = {}

    def register(self, plugin: ToolPlugin) -> None:
        """Register a ToolPlugin. Called at module import time by each tool file."""
        self._plugins[plugin.name] = plugin

    # ------------------------------------------------------------------
    # Query API
    # ------------------------------------------------------------------

    def names(self) -> List[str]:
        return list(self._plugins.keys())

    def categories(self) -> Dict[str, List[str]]:
        cats: Dict[str, List[str]] = defaultdict(list)
        for name, p in self._plugins.items():
            cats[p.category].append(name)
        return dict(cats)

    def get_callables(self) -> Dict[str, Callable[..., Any]]:
        """Return {name: fn} dict for direct tool execution (no LangChain)."""
        return {name: p.fn for name, p in self._plugins.items()}

    def get_langchain_tools(self) -> List[Any]:
        """Return all LangChain @tool objects (for LangGraph integration)."""
        tools = []
        for p in self._plugins.values():
            if p.langchain_tool is not None:
                tools.append(p.langchain_tool)
        return tools

    # ------------------------------------------------------------------
    # Execution group builder (topological sort)
    # ------------------------------------------------------------------

    def get_execution_groups(self, tool_names: List[str]) -> List[List[str]]:
        """
        Sort a list of tool names into sequential groups where tools within
        each group have no inter-dependencies and can run in parallel.

        Example:
          tools = [normalize_reaction, analyze_bond_changes, detect_reaction_type,
                   search_reaction_types, recommend_conditions]
          → Group 0: [normalize_reaction, detect_reaction_type]  (no prerequisites)
          → Group 1: [analyze_bond_changes]                       (needs normalize)
          → Group 2: [search_reaction_types]                      (needs bond_changes)
          → Group 3: [recommend_conditions]                       (needs detect_type)
        """
        # Only consider tools that are actually in the registry
        valid = [t for t in tool_names if t in self._plugins]
        if not valid:
            return []

        # Build prerequisite map restricted to the requested set
        prereqs: Dict[str, List[str]] = {}
        for name in valid:
            plugin = self._plugins[name]
            prereqs[name] = [p for p in plugin.prerequisites if p in valid]

        # Kahn's algorithm for topological sort into levels
        in_degree: Dict[str, int] = {n: len(prereqs[n]) for n in valid}
        dependents: Dict[str, List[str]] = defaultdict(list)
        for name, deps in prereqs.items():
            for dep in deps:
                dependents[dep].append(name)

        groups: List[List[str]] = []
        current_level = [n for n in valid if in_degree[n] == 0]

        while current_level:
            groups.append(sorted(current_level))  # sorted for determinism
            next_level = []
            for name in current_level:
                for dep in dependents[name]:
                    in_degree[dep] -= 1
                    if in_degree[dep] == 0:
                        next_level.append(dep)
            current_level = next_level

        # Any remaining tools with unresolvable deps go to final group
        scheduled = {t for grp in groups for t in grp}
        remaining = [n for n in valid if n not in scheduled]
        if remaining:
            groups.append(sorted(remaining))

        return groups

    # ------------------------------------------------------------------
    # Prompt generation
    # ------------------------------------------------------------------

    def describe_tools(self) -> str:
        """Auto-generate the tool descriptions block for REASON_PROMPT."""
        lines = []
        for cat, names in sorted(self.categories().items()):
            lines.append(f"\n[{cat.upper()}]")
            for name in names:
                p = self._plugins[name]
                prereq_note = ""
                if p.prerequisites:
                    prereq_note = f" (after: {', '.join(p.prerequisites)})"
                lines.append(f"  {name}{prereq_note}")
                lines.append(f"    → {p.description}")
        return "\n".join(lines)

    def __len__(self) -> int:
        return len(self._plugins)

    def __contains__(self, name: str) -> bool:
        return name in self._plugins


# ---------------------------------------------------------------------------
# Global singleton — all tool modules register into this at import time
# ---------------------------------------------------------------------------

REGISTRY = ToolRegistry()


# ---------------------------------------------------------------------------
# Import all tool modules to trigger their REGISTRY.register() calls
# ---------------------------------------------------------------------------

from .chemistry import CHEMISTRY_TOOLS
from .taxonomy import TAXONOMY_TOOLS
from .conditions import CONDITIONS_TOOLS
from .reagent import REAGENT_TOOLS
from .literature import LITERATURE_TOOLS
from .notes import NOTES_TOOLS

for _plugin in CHEMISTRY_TOOLS + TAXONOMY_TOOLS + CONDITIONS_TOOLS + REAGENT_TOOLS + LITERATURE_TOOLS + NOTES_TOOLS:
    REGISTRY.register(_plugin)

# ---------------------------------------------------------------------------
# Convenience: flat list of all registered ToolPlugin objects
# ---------------------------------------------------------------------------

COWORKER_TOOLS = [REGISTRY._plugins[name] for name in REGISTRY.names()]

__all__ = [
    "ToolRegistry",
    "REGISTRY",
    "COWORKER_TOOLS",
]
