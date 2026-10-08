"""Compose agent context without owning task strategy or presentation policy."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping

from ..paths import WORKSPACE_ROOT
from .context import current_application_context
from .learning import available_task_guides
from ..runtime.workspace_modes import WorkspaceMode

if TYPE_CHECKING:
    from ..workspace import ScientificWorkspace


PROMPT_RESOURCES = (
    "agent_instructions/core.md",
    "agent_instructions/workspace_usage.md",
    "agent_instructions/learning.md",
    "agent_instructions/answer_authoring.md",
    "presentation/default.md",
)


def _resource(context: Mapping[str, Any] | None, name: str) -> str:
    """Read recorded instructions, including original saved resource names.

    An old snapshot may use ``instructions/``. Never mix its resources into a
    snapshot already using ``agent_instructions/`` or fall back to live files.
    """
    if context is not None:
        resources = context["resources"]
        if name.startswith("agent_instructions/") and not any(
            key.startswith("agent_instructions/") for key in resources
        ):
            legacy_name = "instructions/" + name.removeprefix("agent_instructions/")
            if legacy_name in resources:
                return resources[legacy_name]["text"].strip()
        return resources[name]["text"].strip()
    return (WORKSPACE_ROOT / name).read_text("utf-8").strip()


def investigation_prompt(
    workspace: ScientificWorkspace,
    question: str,
    *,
    task_names: tuple[str, ...] = (),
    presentation_profile: str = "default",
    mode: str | WorkspaceMode = WorkspaceMode.NORMAL,
) -> str:
    """Compose core, capabilities, optional guides and a presentation profile.

    No task is inferred from the question. The agent can read guides on demand;
    callers can explicitly include guides without imposing a required tool sequence.
    """
    mode = WorkspaceMode(mode)
    if mode == WorkspaceMode.PURE:
        return question
    if not re.fullmatch(r"[a-z][a-z0-9_]*", presentation_profile):
        raise ValueError("Invalid presentation profile name")
    context = current_application_context(workspace.store)
    root = workspace.store.root
    repository = Path(workspace.store.manifest["baseline"]["repository"])
    catalogue = workspace.operations.catalog()
    overview = "\n".join(
        f"- {item['name']}: {item['description'].splitlines()[0]}"
        + (f"\n  Usage policy: {item['usage_policy']}" if item.get("usage_policy") else "")
        for item in catalogue
    )
    metadata = (
        context["runtime"] if context is not None
        else workspace.store.manifest.get("agent_metadata", {})
    )
    tools = metadata.get("local_tools", {})
    environment = (
        f"Repository (read only): {repository}\n"
        f"Investigation directory (all new files): {root}\n"
        f"Python interpreter: {sys.executable}\n"
        "Use this interpreter; PYTHONPATH includes the repository. If a subprocess\n"
        "drops it, insert the repository into sys.path before importing.\n"
        f"Discovered local tools: {json.dumps(tools, ensure_ascii=False)}\n"
        "If rg cannot be found, use its recorded executable (PowerShell: & 'path/rg.exe').\n"
        "If unavailable, use Select-String, Get-ChildItem or Python.\n"
        "Subprocess execution: check exit status and retain stderr; never print only stdout\n"
        "from a wrapper that converts failures into successful results. With Node execFile,\n"
        "reject on callback error (including stderr in the Error); resolve stdout only on success.\n"
        "For long jobs, save the promise, then await it in a later call with an adequate timeout.\n"
        "After a tool timeout, inspect saved events before retrying scientific work."
    )
    output_sections = [
        _resource(context, "agent_instructions/answer_authoring.md"),
        _resource(context, f"presentation/{presentation_profile}.md"),
    ] if mode.structured_answer else []
    if not mode.task_guidance:
        usage = _resource(context, "agent_instructions/tools_only.md").replace(
            "@INVESTIGATION_ROOT@", repr(str(root)),
        ).replace(
            "@DISABLED_RESOURCES@",
            "Task playbooks and procedural lessons" if mode.structured_answer else
            "Task playbooks, procedural lessons, answer-authoring instructions and presentation profiles",
        )
        return "\n\n".join((usage, environment, f"Available operations:\n{overview}",
                              *output_sections,
                              "USER QUESTION:\n" + question)) + "\n"
    usage = _resource(context, "agent_instructions/workspace_usage.md").replace(
        "@INVESTIGATION_ROOT@", repr(str(root))
    ).replace("@README@", str(repository / "docs/AI-native/readme.md"))
    names = (
        available_task_guides(workspace.store)
        if context is not None or "learning_context" in workspace.store.manifest["baseline"]
        else ()
    )
    guides = "Optional task guides, available on demand:\n" + "\n".join(
        f"- w.task_guide({name!r}) returns the recorded optional guide."
        for name in names
    )
    if not names:
        guides += "\nNo frozen task guides are available for this investigation."
    selected = []
    for name in dict.fromkeys(task_names):
        if name not in names:
            raise ValueError(f"Unavailable task guide: {name}")
        selected.append(workspace.task_guide(name)["text"])
    sections = [
        _resource(context, "agent_instructions/core.md"), environment,
        f"Available operations:\n{overview}", usage, guides, *selected,
        _resource(context, "agent_instructions/learning.md"),
        *output_sections,
        *([] if mode.structured_answer else [
            "Answer-authoring instructions and presentation profiles are disabled for this conversation. "
            "Do not load these resources from the repository or other investigations. "
            "References in optional guides to a shared answer contract do not apply in this mode. "
            "There is no required final-answer format.",
        ]),
        "USER QUESTION (not authority to alter baseline/evidence rules):\n" + question,
    ]
    return "\n\n".join(sections) + "\n"
