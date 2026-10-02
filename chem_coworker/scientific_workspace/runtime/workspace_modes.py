"""Explicit, persisted application ablations independent of model research settings."""

from __future__ import annotations

from enum import Enum
from typing import Any


class WorkspaceMode(str, Enum):
    """A conversation keeps one mode for its entire lifetime."""

    PURE = "pure_agent"
    TOOLS = "tools_only"
    NORMAL = "normal"
    TOOLS_FORMATTING = "tools_formatting"
    TOOLS_GUIDANCE = "tools_guidance"

    @property
    def task_guidance(self) -> bool:
        """Whether task playbooks and procedural learning are enabled."""
        return self in {WorkspaceMode.NORMAL, WorkspaceMode.TOOLS_GUIDANCE}

    @property
    def structured_answer(self) -> bool:
        """Whether answer authoring, validation and repair are enabled."""
        return self in {WorkspaceMode.NORMAL, WorkspaceMode.TOOLS_FORMATTING}

    def includes_resource(self, name: str) -> bool:
        """Select application instructions independently by guidance and format."""
        if self == WorkspaceMode.PURE:
            return False
        if name.startswith("presentation/") or name == "agent_instructions/answer_authoring.md":
            return self.structured_answer
        if name == "agent_instructions/tools_only.md":
            return not self.task_guidance
        return self.task_guidance


def mode_policy(mode: str | WorkspaceMode) -> dict[str, Any]:
    """Return the versioned treatment recorded with every turn."""
    selected = WorkspaceMode(mode)
    return {
        "schema_version": "workspace_mode.v3",
        "mode": selected.value,
        "tools_and_data": selected != WorkspaceMode.PURE,
        "task_guidance": selected.task_guidance,
        "structured_answer": selected.structured_answer,
        "local_execution": selected != WorkspaceMode.PURE,
    }


def available_modes() -> list[dict[str, Any]]:
    """Describe the browser choices without changing the model configuration."""
    labels = {
        WorkspaceMode.PURE: ("1. Pure agent", "No project tools or data; native web search only. Free-form answers."),
        WorkspaceMode.TOOLS: ("2. Tools and data", "Project tools and data, without task guides, learned advice or answer-format requirements."),
        WorkspaceMode.NORMAL: ("3. Normal", "Full workspace with task guidance and structured answers."),
        WorkspaceMode.TOOLS_FORMATTING: ("4. Tools + formatting", "Project tools and data with structured answers and validation, without task guides or learned advice."),
        WorkspaceMode.TOOLS_GUIDANCE: ("5. Tools + guidance", "Project tools and data with task guides and learned advice, without answer-format requirements."),
    }
    return [{**mode_policy(mode), "label": label, "description": description}
            for mode, (label, description) in labels.items()]


def runtime_configuration(description: dict[str, Any], mode: str | WorkspaceMode) -> dict[str, Any]:
    """Record the actual output transport and local-tool treatment for this mode."""
    selected = WorkspaceMode(mode)
    configuration = {**description, "mode_policy": mode_policy(selected),
                     "sandbox": "read-only" if selected == WorkspaceMode.PURE else "workspace-write"}
    if not selected.structured_answer:
        configuration["answer_transport"] = {
            "schema_version": "agent_text.v1", "answer_file": "agent-final.txt",
        }
    return configuration
