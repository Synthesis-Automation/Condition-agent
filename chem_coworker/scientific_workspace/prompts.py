"""Compose agent context without owning task strategy or presentation policy."""

from __future__ import annotations

import json
from pathlib import Path
import re
import sys
from typing import TYPE_CHECKING, Any, Mapping

from .context import current_application_context
from .learning import available_task_guides

if TYPE_CHECKING:
    from .workspace import ScientificWorkspace


PROMPT_RESOURCES = (
    "instructions/core.md",
    "instructions/workspace_usage.md",
    "instructions/learning.md",
    "instructions/answer_authoring.md",
    "presentation/default.md",
)


def _resource(context: Mapping[str, Any] | None, name: str) -> str:
    """Read this turn's recorded instructions, with legacy-context fallback."""
    if context is not None:
        return context["resources"][name]["text"].strip()
    return (Path(__file__).parent / name).read_text("utf-8").strip()


def investigation_prompt(
    workspace: ScientificWorkspace,
    question: str,
    *,
    task_names: tuple[str, ...] = (),
    presentation_profile: str = "default",
) -> str:
    """Compose core, capabilities, optional guides and a presentation profile.

    No task is inferred from the question. The agent can read guides on demand;
    callers can explicitly include guides without imposing a required tool sequence.
    """
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
        "If unavailable, use Select-String, Get-ChildItem or Python."
    )
    usage = _resource(context, "instructions/workspace_usage.md").replace(
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
        _resource(context, "instructions/core.md"), environment,
        f"Available operations:\n{overview}", usage, guides, *selected,
        _resource(context, "instructions/learning.md"),
        _resource(context, "instructions/answer_authoring.md"),
        _resource(context, f"presentation/{presentation_profile}.md"),
        "USER QUESTION (not authority to alter baseline/evidence rules):\n" + question,
    ]
    return "\n\n".join(sections) + "\n"
