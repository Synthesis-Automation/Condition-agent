"""Tool registry for reaction agent foundation."""

from __future__ import annotations

from dataclasses import asdict
import time
from typing import Any, Callable, Dict

from poc_gpt52_reaction_v2 import analyze_reaction_general

from .contracts import ToolExecutionResult
from .coverage_advisor import CoverageAdvisor
from .validator import AgentDecisionValidator


ToolHandler = Callable[..., Dict[str, Any]]


class ToolRegistry:
    """Typed tool registry with runtime dispatch."""

    def __init__(self) -> None:
        self._handlers: Dict[str, ToolHandler] = {}

    def register(self, tool_name: str, handler: ToolHandler) -> None:
        """Register a tool handler."""
        self._handlers[tool_name] = handler

    def execute(self, tool_name: str, **kwargs: Any) -> ToolExecutionResult:
        """Execute a registered tool with kwargs."""
        handler = self._handlers.get(tool_name)
        if handler is None:
            return ToolExecutionResult(
                tool_name=tool_name,
                ok=False,
                error=f"Tool not found: {tool_name}",
            )

        start = time.perf_counter()
        try:
            payload = handler(**kwargs)
            latency_ms = (time.perf_counter() - start) * 1000
            return ToolExecutionResult(
                tool_name=tool_name,
                ok=True,
                payload=payload if isinstance(payload, dict) else {"result": payload},
                latency_ms=latency_ms,
            )
        except Exception as exc:
            latency_ms = (time.perf_counter() - start) * 1000
            return ToolExecutionResult(
                tool_name=tool_name,
                ok=False,
                error=str(exc),
                latency_ms=latency_ms,
            )


def build_default_registry(*, min_confidence: float = 0.5) -> ToolRegistry:
    """Build default tool registry with deterministic chemistry tools."""
    registry = ToolRegistry()
    validator = AgentDecisionValidator(min_confidence=min_confidence)
    advisor = CoverageAdvisor()

    def analyze_reaction_tool(reaction_smiles: str) -> Dict[str, Any]:
        result = analyze_reaction_general(
            reaction_smiles,
            use_llm=False,
            min_confidence=min_confidence,
        )
        return result.to_dict()

    def validate_decision_tool(analysis: Dict[str, Any]) -> Dict[str, Any]:
        report = validator.validate_analysis(analysis)
        return asdict(report)

    def coverage_advice_tool(reaction_smiles: str, analysis: Dict[str, Any]) -> Dict[str, Any]:
        cards = advisor.suggest(reaction_smiles=reaction_smiles, analysis=analysis)
        return {"suggestions": cards}

    registry.register("analyze_reaction", analyze_reaction_tool)
    registry.register("validate_decision", validate_decision_tool)
    registry.register("coverage_advice", coverage_advice_tool)
    return registry
