"""
ToolPlugin — self-describing tool unit for the ChemCoworker ToolRegistry.

Every tool in chem_coworker is a ToolPlugin. This makes the system self-describing:
- The ToolRegistry auto-generates tool descriptions for prompts
- Prerequisites enable topological sorting into parallel execution groups
- Adding a new tool = create a ToolPlugin and call REGISTRY.register()
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional


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
    """
    name: str
    category: str
    description: str
    fn: Callable[..., Any]
    prerequisites: List[str] = field(default_factory=list)
    langchain_tool: Optional[Any] = None  # langchain_core.tools.BaseTool

    def __post_init__(self) -> None:
        # If langchain_tool is provided but fn is not explicitly the wrapped function,
        # we store fn separately. fn is always the direct Python callable.
        pass

    def call(self, **kwargs: Any) -> Any:
        """Direct invocation — used by ToolExecutor (no LangChain overhead)."""
        return self.fn(**kwargs)
