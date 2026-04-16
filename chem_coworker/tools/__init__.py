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
from typing import Any, Callable, Dict, List, Optional, Sequence

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

    def filtered_names(
        self,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> List[str]:
        """Return tool names after applying registry-level filters."""
        return self._select_plugin_names(
            llm_exposed_only=llm_exposed_only,
            include_names=include_names,
            exclude_names=exclude_names,
        )

    def _select_plugin_names(
        self,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> List[str]:
        include_list = list(include_names or [])
        include_set = set(include_list)
        exclude_set = set(exclude_names or [])
        selected: List[str] = []
        ordered_names = include_list if include_list else list(self._plugins.keys())
        for name in ordered_names:
            plugin = self._plugins.get(name)
            if plugin is None:
                continue
            if llm_exposed_only and not getattr(plugin, "llm_exposed", True):
                continue
            if include_set and name not in include_set:
                continue
            if name in exclude_set:
                continue
            selected.append(name)
        return selected

    def categories(
        self,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> Dict[str, List[str]]:
        cats: Dict[str, List[str]] = defaultdict(list)
        for name in self._select_plugin_names(
            llm_exposed_only=llm_exposed_only,
            include_names=include_names,
            exclude_names=exclude_names,
        ):
            p = self._plugins[name]
            cats[p.category].append(name)
        return dict(cats)

    def get_callables(self) -> Dict[str, Callable[..., Any]]:
        """Return {name: fn} dict for direct tool execution (no LangChain)."""
        return {name: p.fn for name, p in self._plugins.items()}

    def get_langchain_tools(
        self,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> List[Any]:
        """Return LangChain @tool objects, optionally filtered for LLM exposure."""
        tools = []
        for name in self._select_plugin_names(
            llm_exposed_only=llm_exposed_only,
            include_names=include_names,
            exclude_names=exclude_names,
        ):
            p = self._plugins[name]
            if p.langchain_tool is not None:
                tools.append(p.langchain_tool)
        return tools

    def policy_names(self) -> List[str]:
        return sorted(_TOOL_POLICIES.keys())

    def filtered_names_for_policy(
        self,
        policy_name: str,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> List[str]:
        base_names = list(_TOOL_POLICIES.get(str(policy_name or "").strip(), []))
        if include_names:
            for name in include_names:
                if name not in base_names:
                    base_names.append(name)
        return self._select_plugin_names(
            llm_exposed_only=llm_exposed_only,
            include_names=base_names,
            exclude_names=exclude_names,
        )

    def describe_tools_for_policy(
        self,
        policy_name: str,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> str:
        names = self.filtered_names_for_policy(
            policy_name,
            llm_exposed_only=llm_exposed_only,
            include_names=include_names,
            exclude_names=exclude_names,
        )
        return self.describe_tools(
            llm_exposed_only=llm_exposed_only,
            include_names=names,
            exclude_names=exclude_names,
        )

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

    def describe_tools(
        self,
        *,
        llm_exposed_only: bool = False,
        include_names: Optional[Sequence[str]] = None,
        exclude_names: Optional[Sequence[str]] = None,
    ) -> str:
        """Auto-generate a grouped tool-description block, optionally filtered."""
        lines = []
        categories = self.categories(
            llm_exposed_only=llm_exposed_only,
            include_names=include_names,
            exclude_names=exclude_names,
        )
        for cat, names in sorted(categories.items()):
            lines.append(f"\n[{cat.upper()}]")
            for name in sorted(names):
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
from .retrosynthesis import RETROSYNTHESIS_TOOLS
from .name_resolver import NAME_RESOLVER_TOOLS
from .reaction_eval import EVAL_TOOLS
from .molecular_features import MOLECULAR_FEATURE_TOOLS
from .forward_synthesis import FORWARD_SYNTHESIS_TOOLS
from .composite import COMPOSITE_TOOLS

for _plugin in (
    CHEMISTRY_TOOLS + TAXONOMY_TOOLS + CONDITIONS_TOOLS
    + REAGENT_TOOLS + LITERATURE_TOOLS + NOTES_TOOLS
    + RETROSYNTHESIS_TOOLS + NAME_RESOLVER_TOOLS + EVAL_TOOLS
    + MOLECULAR_FEATURE_TOOLS + FORWARD_SYNTHESIS_TOOLS + COMPOSITE_TOOLS
):
    REGISTRY.register(_plugin)


# ---------------------------------------------------------------------------
# Default LLM exposure profile (facade-first)
# ---------------------------------------------------------------------------
#
# Keep all tools registered and executable, but hide most legacy primitive tools
# from the default LLM/public view. Workflow-specific allowlists already curate
# the actual task-time tool surface; this primarily reduces /tools noise and
# accidental exposure in generic contexts.
_HIDE_FROM_LLM_BY_DEFAULT = {
    # Primitive support tools replaced by composite facades
    "resolve_to_smiles",
    "smiles_to_info",
    "lookup_reagent",
    "list_reagents_by_role",
    "search_literature",
    "fetch_webpage",
    "read_literature_source",
    "read_notes",
    "search_notes",
    "list_notes",
    "search_knowledge",
    "read_knowledge",
    # Primitive reaction analysis / validation replaced by facades
    "detect_reaction_type",
    "inspect_functional_groups",
    "get_molecular_descriptors",
    "recommend_conditions",
    "evaluate_reaction",
    "check_retro_consistency",
    # Primitive retrosynthesis steps replaced by retrosynthesis_step facade
    "inspect_target",
    "identify_retrons",
    "generate_disconnections",
    "find_retro_precedent",
    "search_hte_precedent",
    # Primitive forward synthesis steps replaced by forward_synthesis_step facade
    "inspect_reactants",
    "identify_reactions",
    "generate_products",
    "rank_products",
    "find_forward_precedent",
    "search_reactant_precedent",
    "recommend_forward_conditions",
}

for _name in _HIDE_FROM_LLM_BY_DEFAULT:
    if _name in REGISTRY._plugins:
        REGISTRY._plugins[_name].llm_exposed = False


# ---------------------------------------------------------------------------
# Named tool policies
# ---------------------------------------------------------------------------

_TOOL_POLICIES: Dict[str, List[str]] = {
    "general_chemistry": [
        "analyze_reaction",
        "retrosynthesis_step",
        "forward_synthesis_step",
        "plan_route",
        "plan_forward_route",
        "evaluate_synthesis_proposal",
        "validate_synthesis_proposal",
        "featurize_molecule",
        "assess_snar_feasibility",
        "resolve_chemical",
        "reagent_assistant",
        "recommend_reaction_conditions",
    ],
    "retro_specialist": [
        "retrosynthesis_step",
        "plan_route",
        "apply_hte_templates",
        "search_by_product_similarity",
        "evaluate_synthesis_proposal",
        "validate_synthesis_proposal",
        "featurize_molecule",
        "assess_snar_feasibility",
        "resolve_chemical",
        "reagent_assistant",
        "recommend_reaction_conditions",
    ],
    "forward_specialist": [
        "forward_synthesis_step",
        "plan_forward_route",
        "evaluate_synthesis_proposal",
        "validate_synthesis_proposal",
        "featurize_molecule",
        "assess_snar_feasibility",
        "resolve_chemical",
        "reagent_assistant",
        "recommend_reaction_conditions",
    ],
    "conditions_specialist": [
        "build_reaction_context",
        "detect_reaction_type",
        "inspect_functional_groups",
        "get_literature_condition_evidence",
        "get_motif_condition_evidence",
        "get_similarity_condition_evidence",
        "get_rule_condition_evidence",
        "compose_condition_candidates",
        "analyze_reaction",
        "recommend_reaction_conditions",
        "reagent_assistant",
        "featurize_molecule",
        "resolve_chemical",
    ],
    "literature_curator": [
        "search_notes",
        "read_notes",
        "search_knowledge",
        "read_knowledge",
    ],
    "internal_debug": REGISTRY.names(),
}

# ---------------------------------------------------------------------------
# Convenience: flat list of all registered ToolPlugin objects
# ---------------------------------------------------------------------------

COWORKER_TOOLS = [REGISTRY._plugins[name] for name in REGISTRY.names()]

__all__ = [
    "ToolRegistry",
    "REGISTRY",
    "COWORKER_TOOLS",
    "_TOOL_POLICIES",
]
