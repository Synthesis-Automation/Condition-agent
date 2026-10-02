"""Prompt boundaries keep task advice and display policy outside orchestration."""

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from chem_coworker.scientific_workspace import prompts


@pytest.fixture
def prompt_context(tmp_path, monkeypatch):
    directory = Path(prompts.__file__).parent
    resources = {
        name: {"text": (directory / name).read_text("utf-8")}
        for name in prompts.PROMPT_RESOURCES
    }
    context = {
        "resources": resources,
        "runtime": {"local_tools": {"rg": {"path": "recorded/rg.exe"}}},
    }
    guides = {"new_task": {"text": "An optional new-task strategy."}}
    workspace = SimpleNamespace(
        store=SimpleNamespace(root=tmp_path, manifest={
            "baseline": {"repository": str(tmp_path)},
            "agent_metadata": {"local_tools": {"rg": {"path": "old/rg.exe"}}},
        }),
        operations=SimpleNamespace(catalog=lambda: [{
            "name": "new_operation", "description": "A new scientific capability.",
            "usage_policy": "Respect this capability's evidence limitations.",
        }]),
        task_guide=lambda name: deepcopy(guides[name]),
    )
    monkeypatch.setattr(prompts, "current_application_context", lambda _: context)
    monkeypatch.setattr(prompts, "available_task_guides", lambda _: tuple(guides))
    return workspace, context, guides


def test_conversation_reexports_single_prompt_implementation():
    from chem_coworker.scientific_workspace.conversation import investigation_prompt

    assert investigation_prompt is prompts.investigation_prompt


def test_core_discovers_new_capabilities_and_guides_without_task_routing(prompt_context):
    workspace, _, guides = prompt_context
    text = prompts.investigation_prompt(workspace, "Plan a synthesis and conditions")
    assert "new_operation: A new scientific capability." in text
    assert "Usage policy: Respect this capability's evidence limitations." in text
    assert "w.task_guide('new_task')" in text
    assert guides["new_task"]["text"] not in text
    assert "recorded/rg.exe" in text and "old/rg.exe" not in text
    core = text.split("Repository (read only):")[0]
    for task_api in ("disconnect_target", "recommend_conditions", "finalize_answer"):
        assert task_api not in core
    assert "choose your own" in core.lower()


def test_recorded_empty_runtime_does_not_revive_stale_tool_metadata(prompt_context):
    workspace, context, _ = prompt_context
    context["runtime"] = {}
    text = prompts.investigation_prompt(workspace, "Investigate")
    assert "old/rg.exe" not in text


def test_explicit_guidance_is_optional_deduplicated_and_frozen(prompt_context, monkeypatch):
    workspace, _, guides = prompt_context
    # An application context uses snapshots, even if current files cannot be read.
    monkeypatch.setattr(Path, "read_text", lambda *_, **__: pytest.fail("Unrecorded resource read"))
    text = prompts.investigation_prompt(
        workspace, "Investigate", task_names=("new_task", "new_task"),
    )
    assert text.count(guides["new_task"]["text"]) == 1
    assert "optional menus" in text
    with pytest.raises(ValueError, match="Unavailable task guide"):
        prompts.investigation_prompt(workspace, "Investigate", task_names=("unknown",))


def test_presentation_changes_independently_of_core_or_task_guide(prompt_context):
    workspace, context, _ = prompt_context
    first = prompts.investigation_prompt(workspace, "Investigate")
    context["resources"]["presentation/default.md"]["text"] = "Use a short text answer."
    second = prompts.investigation_prompt(workspace, "Investigate")
    assert second.split("Use a short text answer.")[0] == first.split(
        "Default answer presentation"
    )[0]
    assert "scientific_answer.v2" in second
    assert "matching\nnonempty inspect_step_precedents inspection" in second


def test_new_context_never_silently_uses_missing_current_resource(prompt_context):
    workspace, context, _ = prompt_context
    del context["resources"]["instructions/core.md"]
    with pytest.raises(KeyError, match="instructions/core.md"):
        prompts.investigation_prompt(workspace, "Investigate")
    with pytest.raises(ValueError, match="presentation profile"):
        prompts.investigation_prompt(workspace, "Investigate", presentation_profile="../core")


def test_legacy_investigation_without_guidance_can_still_compose(tmp_path, monkeypatch):
    workspace = SimpleNamespace(
        store=SimpleNamespace(root=tmp_path, manifest={"baseline": {"repository": str(tmp_path)}}),
        operations=SimpleNamespace(catalog=lambda: []),
    )
    monkeypatch.setattr(prompts, "current_application_context", lambda _: None)
    monkeypatch.setattr(prompts, "available_task_guides", lambda _: pytest.fail("No legacy guides"))
    text = prompts.investigation_prompt(workspace, "A general question")
    assert "No frozen task guides" in text
    assert "A general question" in text
    assert "w.run_python" in text
