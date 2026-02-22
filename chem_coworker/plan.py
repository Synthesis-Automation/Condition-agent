"""
Core data types for ChemCoworker's native tool-calling pipeline.

The classic text-plan path (PlanParser + JSON plan) was removed in Phase 5.
The native tool-calling path uses LangChain's bind_tools() and ToolMessage API
directly; these data classes are retained for ToolExecutor's parallel execution
engine and for the plan_callback approval gate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List


class PlanRejected(Exception):
    """Raised by a plan_callback to cancel tool execution (A2 — plan approval)."""


@dataclass
class ToolCall:
    """A single tool call specification (name + arguments)."""
    name: str
    args: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        args_str = ", ".join(f"{k}={v!r}" for k, v in self.args.items())
        return f"{self.name}({args_str})"


@dataclass
class ExecutionPlan:
    """
    Minimal execution plan used by the plan_callback approval gate.

    In the native tool-calling path the LLM drives tool selection directly;
    ExecutionPlan is populated from the native loop's final tool list so the
    plan_callback can still inspect and approve/reject tool usage.

    Attributes:
        hypothesis: Named reaction / chemistry concept.
        confidence: 0.0–1.0 reflecting LLM certainty.
        groups: Sequential groups of ToolCalls (kept for plan_callback display).
        rationale: Brief explanation of the approach.
        raw_plan_text: Unused in native mode; kept for backward compatibility.
    """
    hypothesis: str = ""
    confidence: float = 0.5
    groups: List[List[ToolCall]] = field(default_factory=list)
    rationale: str = ""
    raw_plan_text: str = ""

    @property
    def all_tool_calls(self) -> List[ToolCall]:
        return [tc for grp in self.groups for tc in grp]

    @property
    def tool_names(self) -> List[str]:
        return [tc.name for tc in self.all_tool_calls]

    @property
    def is_empty(self) -> bool:
        """True if no tools are planned (pure LLM reasoning task)."""
        return len(self.all_tool_calls) == 0

    def __repr__(self) -> str:
        groups_repr = " → ".join(
            f"[{', '.join(tc.name for tc in grp)}]" for grp in self.groups
        )
        return (
            f"ExecutionPlan(hypothesis={self.hypothesis!r}, "
            f"confidence={self.confidence:.2f}, groups={groups_repr})"
        )
