"""Explicit, persisted application ablations independent of model research settings."""

from __future__ import annotations

from enum import Enum
from typing import Any


class WorkspaceMode(str, Enum):
    """A conversation keeps one mode for its entire lifetime."""

    PURE = "pure_agent"
    TOOLS = "tools_only"
    NORMAL = "normal"


def mode_policy(mode: str | WorkspaceMode) -> dict[str, Any]:
    """Return the versioned treatment recorded with every turn."""
    selected = WorkspaceMode(mode)
    return {
        "schema_version": "workspace_mode.v2",
        "mode": selected.value,
        "tools_and_data": selected != WorkspaceMode.PURE,
        "task_guidance": selected == WorkspaceMode.NORMAL,
        "structured_answer": selected == WorkspaceMode.NORMAL,
        "local_execution": selected != WorkspaceMode.PURE,
    }


def available_modes() -> list[dict[str, Any]]:
    """Describe the three browser choices without changing the model configuration."""
    labels = {
        WorkspaceMode.PURE: ("1. Pure agent", "No project tools or data; native web search only. Free-form answers."),
        WorkspaceMode.TOOLS: ("2. Tools and data", "Project tools and data, without task guides, learned advice or answer-format requirements."),
        WorkspaceMode.NORMAL: ("3. Normal", "Full workspace with task guidance and structured answers."),
    }
    return [{**mode_policy(mode), "label": label, "description": description}
            for mode, (label, description) in labels.items()]


def runtime_configuration(description: dict[str, Any], mode: str | WorkspaceMode) -> dict[str, Any]:
    """Record the actual output transport and local-tool treatment for this mode."""
    selected = WorkspaceMode(mode)
    if selected == WorkspaceMode.NORMAL:
        return dict(description)
    return {**description, "mode_policy": mode_policy(selected),
            "answer_transport": {"schema_version": "agent_text.v1", "answer_file": "agent-final.txt"},
            "sandbox": "read-only" if selected == WorkspaceMode.PURE else "workspace-write"}
