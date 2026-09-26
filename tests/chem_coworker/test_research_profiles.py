"""Execution requests, CLI precedence, and observed runtime activity regressions."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from threading import Event

import pytest

from chem_coworker.scientific_workspace.agent_runtime import CodexRuntime
from chem_coworker.scientific_workspace.research_profiles import resolve_research_profile


@pytest.fixture
def local_runtime(monkeypatch: pytest.MonkeyPatch):
    """Avoid discovering or authenticating a real model runtime in unit tests."""
    monkeypatch.setattr(
        "chem_coworker.scientific_workspace.agent_runtime.find_codex",
        lambda executable: (executable or "codex-test", "test-version"),
    )
    return CodexRuntime


@pytest.mark.parametrize("profile,reasoning,search,timeout", [
    ("inherit", None, None, 900),
    ("quick", "medium", "cached", 300),
    ("research", "high", "live", 1800),
])
def test_profiles_leave_model_selection_unchanged(local_runtime, profile, reasoning, search, timeout) -> None:
    runtime = local_runtime(profile=profile)
    assert runtime.model is None
    assert runtime.settings.reasoning_effort == reasoning
    assert runtime.settings.web_search == search
    assert runtime.timeout_seconds == timeout
    description = runtime.describe()
    assert description["model_request"] == {"value": None, "source": "inherited"}
    assert description["effective_configuration_status"] == "unconfirmed"
    assert set(description["capabilities"].values()) == {"not_checked"}


def test_explicit_overrides_win_without_silent_model_substitution(local_runtime, tmp_path: Path) -> None:
    runtime = local_runtime(
        profile="research", model="chosen-provider-model", reasoning_effort="low",
        web_search="disabled", timeout_seconds=123.5,
    )
    assert runtime.settings.to_dict()["sources"] == {
        "reasoning_effort": "override", "web_search": "override", "timeout_seconds": "override",
    }
    assert runtime.timeout_seconds == 123.5
    command = runtime.command(tmp_path, tmp_path, None)
    assert command[command.index("--model") + 1] == "chosen-provider-model"
    assert 'model_reasoning_effort="low"' in command
    assert 'web_search="disabled"' in command
    assert "high" not in command


@pytest.mark.parametrize("thread", [None, "specific-thread"])
def test_config_overrides_apply_to_new_and_resumed_turns(local_runtime, tmp_path: Path, thread) -> None:
    runtime = local_runtime(profile="research")
    command = runtime.command(tmp_path, tmp_path / "turn", thread)
    assert command[:4] == ["codex-test", "exec", "--sandbox", "workspace-write"]
    assert 'approval_policy="never"' in command
    for override in ('model_reasoning_effort="high"', 'web_search="live"'):
        index = command.index(override)
        assert command[index - 1] == "-c"
        if thread:
            assert index > command.index("resume")
    if thread:
        assert command[command.index("resume") + 1] == thread
        assert command.index("--cd") < command.index("resume")
    assert "--model" not in command
    assert "--last" not in command
    assert command[-1] == "-"


def test_default_constructor_retains_inherited_codex_settings(local_runtime, tmp_path: Path) -> None:
    runtime = local_runtime()
    command = runtime.command(tmp_path, tmp_path, None)
    assert runtime.timeout_seconds == 900
    assert not any(value.startswith(("model_reasoning_effort=", "web_search=")) for value in command)
    assert runtime.settings.to_dict()["sources"]["reasoning_effort"] == "inherited"


@pytest.mark.parametrize("kwargs", [
    {"profile": "unknown"}, {"reasoning_effort": "invented"}, {"web_search": "on"},
    {"timeout_seconds": float("nan")}, {"timeout_seconds": float("inf")},
    {"timeout_seconds": True}, {"timeout_seconds": 0}, {"timeout_seconds": 7201},
])
def test_invalid_profile_settings_fail_before_runtime_discovery(kwargs) -> None:
    with pytest.raises(ValueError):
        resolve_research_profile(**kwargs)


def test_runtime_records_requested_configuration_and_unique_observed_events(
    local_runtime, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    runtime = local_runtime(profile="research", model="requested-model")
    events = [
        {"type": "thread.started", "thread_id": "observed-thread"},
        {"type": "item.started", "item": {"id": "search-1", "type": "web_search", "status": "in_progress"}},
        {"type": "item.completed", "item": {"id": "search-1", "type": "web_search", "status": "completed"}},
        {"type": "item.completed", "item": {"id": "search-1", "type": "web_search", "status": "completed"}},
        {"type": "item.completed", "item": {"id": "search-2", "type": "web_search", "status": "failed"}},
        {"type": "item.completed", "item": {"id": "code-1", "type": "command_execution", "status": "completed"}},
        {"type": "turn.completed", "usage": {"input_tokens": 7}},
    ]
    script = tmp_path / "runtime_double.py"
    script.write_text(
        "import json, pathlib, sys\n"
        "pathlib.Path('agent-final.json').write_text(json.dumps({'answer_markdown': 'test'}))\n"
        f"events = json.loads({json.dumps(json.dumps(events))})\n"
        "sys.stdout.write('\\n'.join(json.dumps(event) for event in events))\n",
        "utf-8",
    )
    monkeypatch.setattr(runtime, "command", lambda *_: [sys.executable, str(script)])
    progress = []
    result = runtime.run(
        prompt="Test", workspace=tmp_path, turn_directory=tmp_path,
        thread_id="prior-thread", cancel=Event(), on_event=progress.append,
    )
    request = json.loads((tmp_path / "runtime-request.json").read_text("utf-8"))
    observed = json.loads((tmp_path / "runtime-observations.json").read_text("utf-8"))
    assert request["resumed_thread_id"] == "prior-thread"
    assert request["configuration"]["model_request"]["value"] == "requested-model"
    assert observed["tool_events"]["web_search"] == {
        "observed": 2, "completed_events": 2, "reported_failed": 1,
    }
    assert observed["tool_events"]["command_execution"]["observed"] == 1
    assert observed["tool_events"]["mcp_tool_call"]["observed"] == 0
    assert observed["turn_completed_event"] is True
    assert result.thread_id == "observed-thread"
    assert result.usage == {"input_tokens": 7}
    assert len(progress) == len(events)
    assert "effective_model" not in observed


def test_failed_runtime_retains_request_and_observation_metadata(
    local_runtime, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    runtime = local_runtime(profile="quick")
    script = tmp_path / "failed_runtime.py"
    script.write_text("import sys\nprint('{\"type\":\"turn.failed\"}')\nsys.exit(3)\n", "utf-8")
    monkeypatch.setattr(runtime, "command", lambda *_: [sys.executable, str(script)])
    with pytest.raises(RuntimeError, match="exit 3"):
        runtime.run(prompt="Test", workspace=tmp_path, turn_directory=tmp_path,
                    thread_id=None, cancel=Event(), on_event=lambda _: None)
    observed = json.loads((tmp_path / "runtime-observations.json").read_text("utf-8"))
    assert observed["turn_failed_event"] is True
    assert observed["turn_completed_event"] is False
    assert observed["process_exit_code"] == 3
    assert (tmp_path / "runtime-request.json").is_file()


@pytest.mark.parametrize("flags,expected", [
    ([], ("research", "high", "live", 1800, None)),
    (["--agent-profile", "inherit"], ("inherit", None, None, 900, None)),
    (["--agent-profile", "quick", "--agent-reasoning-effort", "high",
      "--agent-web-search", "disabled", "--agent-timeout", "42", "--agent-model", "chosen"],
     ("quick", "high", "disabled", 42, "chosen")),
])
def test_web_cli_profile_defaults_and_override_precedence(
    local_runtime, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, flags, expected,
) -> None:
    from app.web_api import __main__ as web_main

    captured = {}

    class Service:
        def __init__(self, root, *, runtime, artifacts):
            captured["runtime"] = runtime

        def close(self):
            captured["closed"] = True

    artifacts = tmp_path / "artifacts.json"
    artifacts.write_text("{}", "utf-8")
    monkeypatch.setattr("chem_coworker.scientific_workspace.conversation.ConversationService", Service)
    monkeypatch.setattr(web_main, "LocalRecommendationRuntime", lambda _: object())
    monkeypatch.setattr(web_main, "create_app", lambda **_: object())
    monkeypatch.setattr(web_main.uvicorn, "run", lambda *_, **__: None)
    monkeypatch.setattr(sys, "argv", ["app.web_api", "--scientific-chat", "--chat-artifacts", str(artifacts), *flags])
    web_main.main()
    runtime = captured["runtime"]
    assert (runtime.settings.profile, runtime.settings.reasoning_effort,
            runtime.settings.web_search, runtime.timeout_seconds, runtime.model) == expected
    assert captured["closed"] is True
