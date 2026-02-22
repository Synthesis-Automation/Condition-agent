"""
ChemResponse — the output dataclass for ChemCoworker.

Captures everything from a single run:
  - The user's query and its classified task type
  - The LLM's reasoning and hypothesis (Phase 1)
  - All tool results (Phase 2)
  - The final synthesized answer (Phase 3)
  - Performance metadata (timing, LLM call count, model used)
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class ChemResponse:
    """
    Complete output of a ChemCoworker.run() call.

    Key fields:
        query: The original user query.
        task_type: Detected task type (analyze/predict/explain/lookup/compare/troubleshoot).
        hypothesis: LLM's named reaction or chemistry concept from Phase 1.
        plan_rationale: Why the LLM chose the tool plan it did.
        answer: The final synthesized answer from Phase 3 (human-readable).
        tools_called: Ordered list of tool names that were actually executed.
        tool_results: Raw payload from each tool, keyed by tool name.
        structured: Machine-readable key outputs extracted from tool results
                    (e.g. reaction_type, conditions, descriptors).
        confidence: LLM self-reported confidence in the hypothesis/answer.
        warnings: Any issues encountered (tool failures, parse errors, etc.).
        model: LLM model name used.
        elapsed_s: Wall-clock time for the full run.
        llm_calls: Number of LLM calls made (typically 2).
    """
    # Input
    query: str = ""
    task_type: str = "general"

    # Phase 1 — Reasoning
    hypothesis: str = ""
    plan_rationale: str = ""

    # Phase 3 — Answer
    answer: str = ""

    # Execution record
    tools_called: List[str] = field(default_factory=list)
    tool_results: Dict[str, Any] = field(default_factory=dict)

    # Extracted machine-readable outputs
    structured: Dict[str, Any] = field(default_factory=dict)

    # Metadata
    confidence: float = 0.0
    warnings: List[str] = field(default_factory=list)
    model: str = ""
    provider: str = ""
    elapsed_s: float = 0.0
    llm_calls: int = 0
    compacted: bool = False   # True if conversation history was compacted this turn (A4)
    streamed: bool = False    # True if answer was written to stdout token-by-token

    # ------------------------------------------------------------------
    # Convenience methods
    # ------------------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a plain dictionary."""
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Serialize to JSON string."""
        return json.dumps(self.to_dict(), indent=indent, default=str)

    @property
    def success(self) -> bool:
        """True if an answer was produced without critical failures."""
        return bool(self.answer) and "critical failure" not in self.answer.lower()

    @property
    def reaction_type(self) -> Optional[str]:
        """Shortcut to reaction_type from structured output."""
        return self.structured.get("reaction_type")

    @property
    def conditions(self) -> List[Dict[str, Any]]:
        """Shortcut to recommended conditions from structured output."""
        return self.structured.get("conditions", [])

    def summary(self) -> str:
        """One-line summary for logging/display."""
        parts = [f"task={self.task_type}"]
        if self.hypothesis:
            h = self.hypothesis[:60] + "..." if len(self.hypothesis) > 60 else self.hypothesis
            parts.append(f"hypothesis={h!r}")
        if self.compacted:
            parts.append("compacted=True")
        if self.tools_called:
            parts.append(f"tools={self.tools_called}")
        parts.append(f"conf={self.confidence:.2f}")
        parts.append(f"t={self.elapsed_s:.1f}s")
        parts.append(f"llm_calls={self.llm_calls}")
        return " | ".join(parts)
