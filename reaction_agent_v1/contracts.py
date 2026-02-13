"""Core contracts for the reaction agent foundation."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional


def utc_now_iso() -> str:
    """Return UTC timestamp in ISO 8601 format."""
    return datetime.now(timezone.utc).isoformat()


@dataclass
class PlanStep:
    """A single planner decision."""

    action: str
    reason: str
    tool_name: Optional[str] = None


@dataclass
class ToolExecutionResult:
    """Result of a tool invocation."""

    tool_name: str
    ok: bool
    payload: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    latency_ms: float = 0.0


@dataclass
class TraceEvent:
    """One trace event in a session run."""

    step_index: int
    timestamp: str
    action: str
    tool_name: Optional[str]
    status: str
    message: str


@dataclass
class ValidationReport:
    """Validator output for agent-level decision checks."""

    passed: bool
    issues: List[str] = field(default_factory=list)
    final_decision: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CoverageSuggestionCard:
    """One suggestion for taxonomy or tool coverage expansion."""

    suggestion_id: str
    kind: str
    priority: int
    title: str
    rationale: str
    proposed_changes: List[str] = field(default_factory=list)
    evidence: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SessionState:
    """State container for one session lane."""

    session_id: str
    reaction_smiles: str
    status: str = "running"
    context: Dict[str, Any] = field(default_factory=dict)
    trace: List[TraceEvent] = field(default_factory=list)


@dataclass
class AgentRunResult:
    """Top-level output from a gateway run."""

    session_id: str
    reaction_smiles: str
    status: str
    final_decision: Dict[str, Any]
    evidence: Dict[str, Any]
    analysis: Dict[str, Any]
    validation: Dict[str, Any]
    tool_artifacts: Dict[str, Any]
    coverage_suggestions: List[Dict[str, Any]]
    trace: List[TraceEvent]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-friendly dictionary."""
        return asdict(self)
