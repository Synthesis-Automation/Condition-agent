"""
ToolPlugin — self-describing tool unit for the ChemCoworker ToolRegistry.

Every tool in chem_coworker is a ToolPlugin. This makes the system self-describing:
- The ToolRegistry auto-generates tool descriptions for prompts
- Prerequisites enable topological sorting into parallel execution groups
- Adding a new tool = create a ToolPlugin and call REGISTRY.register()
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional


# Validator signature: fn(result_dict) -> Optional[warning_string]
# Return None = pass, return str = warning surfaced in ChemResponse
ValidatorFn = Callable[[Dict[str, Any]], Optional[str]]


@dataclass
class ToolPlugin:
    """A self-describing, registerable tool unit.

    Attributes:
        name: Unique tool name (matches the @tool function name exactly)
        category: Logical grouping ('chemistry', 'taxonomy', 'conditions', 'reagent', ...)
        description: One-line description shown to the LLM in the REASON_PROMPT
        prerequisites: Names of tools that must complete before this one can run.
            Used by ToolRegistry.get_execution_groups() to build parallel batches.
        fn: The raw callable (Python function, not LangChain tool) for direct execution.
            This is what the ToolExecutor calls — no LangChain overhead.
        langchain_tool: Optional LangChain @tool object, for direct LangGraph use.
        provides: Keys this tool writes into the shared result context.
            Examples: ["resolved_smiles"], ["reaction_type"], ["conditions"]
            Used by executor for data-contract resolution (replaces _SMILES_PLACEHOLDER_RE).
        requires: Keys from prior results this tool needs (informational; executor
            logs a warning if missing but does not abort).
        validators: Post-execution validators run by ToolExecutor after the tool returns.
            Each fn(result_dict) -> Optional[str]. Non-None return is a warning
            appended to result["_warnings"] and surfaced in ChemResponse.
            Replaces hardcoded tool-name checks in agent._check_hypothesis().
    """
    name: str
    category: str
    description: str
    fn: Callable[..., Any]
    prerequisites: List[str] = field(default_factory=list)
    langchain_tool: Optional[Any] = None  # langchain_core.tools.BaseTool
    # --- Data-contract metadata (Phase 1) ---
    provides: List[str] = field(default_factory=list)
    requires: List[str] = field(default_factory=list)
    validators: List[ValidatorFn] = field(default_factory=list)

    def __post_init__(self) -> None:
        # If langchain_tool is provided but fn is not explicitly the wrapped function,
        # we store fn separately. fn is always the direct Python callable.
        pass

    def call(self, **kwargs: Any) -> Any:
        """Direct invocation — used by ToolExecutor (no LangChain overhead)."""
        return self.fn(**kwargs)

    def to_langchain_tool(self) -> Any:
        """Create a LangChain StructuredTool from this ToolPlugin.

        Schema is inferred automatically from the underlying function's type
        annotations and docstring. Used by the native tool-calling loop in
        ChemCoworker when native_tools=True.
        """
        from langchain_core.tools import StructuredTool
        return StructuredTool.from_function(
            func=self.fn,
            name=self.name,
            description=self.description,
        )
