"""Ablation boundaries, free-form transport and persistent browser mode selection."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from threading import Event
from time import monotonic, sleep
from types import SimpleNamespace

from fastapi import FastAPI
from fastapi.testclient import TestClient
import pytest

from app.web_api.scientific_chat import create_scientific_router
from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.agent_context.context import current_application_context
from chem_coworker.scientific_workspace.agent_context.learning import (
    DEVELOPMENT_PARTITION, build_learning_context, disabled_learning_context,
)
from chem_coworker.scientific_workspace.answers.answer_handoff import load_text_answer, prepare_answer_handoff
from chem_coworker.scientific_workspace.runtime.agent_runtime import AgentResult, CodexRuntime
from chem_coworker.scientific_workspace.runtime.conversation import ConversationService
from chem_coworker.scientific_workspace.runtime.research_profiles import resolve_research_profile
from chem_coworker.scientific_workspace.runtime.runtime_environment import runtime_environment
from chem_coworker.scientific_workspace.runtime.workspace_modes import WorkspaceMode, mode_policy


class FreeRuntime:
    """Capture prompts and simulate unrestricted final text, without a provider."""

    def __init__(self) -> None:
        self.calls: list[dict] = []

    def describe(self) -> dict:
        return {"runtime": "test", "model": "fixed-model"}

    def run(self, **kwargs) -> AgentResult:
        self.calls.append(kwargs)
        return AgentResult({"answer_markdown": "Unstructured text, no JSON required.\n"}, "thread", {"input_tokens": 3})


def finish(service: ConversationService, identity: str) -> dict:
    deadline = monotonic() + 15
    while monotonic() < deadline:
        if service._active is None:
            return service.get(identity)["turns"][-1]
        sleep(.01)
    raise AssertionError("Worker did not finish")


@pytest.fixture
def service(tmp_path, monkeypatch):
    repository = tmp_path / "repository"
    repository.mkdir()

    def baseline(repo, artifacts, *, include_guidance=True):
        result = {"repository": str(repo), "code_files": {}, "environment": {}, "artifacts": {},
                  "evaluation_partition": DEVELOPMENT_PARTITION}
        result["learning_context"] = (build_learning_context(result) if include_guidance
                                      else disabled_learning_context())
        return result

    monkeypatch.setattr("chem_coworker.scientific_workspace.core.baseline.capture_baseline", baseline)
    monkeypatch.setattr("chem_coworker.scientific_workspace.core.baseline.environment_versions", lambda: {})
    instance = ConversationService(tmp_path / "chats", runtime=FreeRuntime(), repository=repository)
    yield instance
    instance.close()


@pytest.mark.parametrize("mode", [WorkspaceMode.PURE, WorkspaceMode.TOOLS, WorkspaceMode.TOOLS_GUIDANCE])
def test_free_answers_and_mode_survive_followup_and_restart(service, mode):
    first = service.submit("What can you infer?", mode=mode)
    identity = first["conversation_id"]
    turn = finish(service, identity)
    assert turn["status"] == "completed", turn.get("error")
    assert turn["answer"]["schema_version"] == "agent_text.v1"
    assert turn["answer"]["answer_markdown"] == "Unstructured text, no JSON required.\n"
    assert turn["answer"]["evidence_status"] == "not_validated"
    assert "repair_attempts" not in turn
    assert turn["workspace_mode"] == mode.value
    assert service.runtime.calls[0]["mode"] == mode
    assert turn["runtime_requested"]["answer_transport"]["answer_file"] == "agent-final.txt"
    with pytest.raises(ValueError, match="fixed per conversation"):
        service.submit("Change treatment", identity, mode="normal")
    reopened = ConversationService(service.root, runtime=service.runtime, repository=service.repository)
    try:
        reopened.submit("Follow up", identity)
        second = finish(reopened, identity)
        assert second["status"] == "completed", second.get("error")
        assert second["thread_resume"] == "resumed"
        assert service.runtime.calls[-1]["thread_id"] == "thread"
        assert reopened.list_conversations()[0]["workspace_mode"] == mode.value
    finally:
        reopened.close()


def test_pure_has_no_scientific_context_or_baseline(service, monkeypatch):
    monkeypatch.setattr("chem_coworker.scientific_workspace.core.baseline.capture_baseline",
                        lambda *a, **k: pytest.fail("Pure mode must not access datasets"))
    identity = service.submit("Exactly this question", mode="pure_agent")["conversation_id"]
    turn = finish(service, identity)
    assert turn["status"] == "completed"
    assert service.runtime.calls[0]["prompt"] == "Exactly this question"
    assert not (service.root / identity / "investigation.json").exists()
    assert "application_context_ref" not in turn
    assert "scientific_identity" not in turn


def test_tools_prompt_and_saved_context_exclude_suggestion_and_output_layers(service):
    identity = service.submit("Investigate", mode="tools_only")["conversation_id"]
    turn = finish(service, identity)
    assert turn["status"] == "completed", turn.get("error")
    prompt = service.runtime.calls[0]["prompt"]
    assert "analyze_reaction" in prompt and "w.run_summary(operation_name, parameters)" in prompt
    for forbidden in ("w.task_guide(", "w.recall_lessons(", "scientific_answer.v2",
                      "ANSWER FILE HANDOFF", "Default answer presentation", "optional menus"):
        assert forbidden not in prompt
    workspace = ScientificWorkspace(service.root / identity)
    context = current_application_context(workspace.store)
    assert set(context["resources"]) == {"agent_instructions/tools_only.md"}
    assert context["learning_context"] == disabled_learning_context()
    assert workspace.recall_lessons("conditions")["lessons"] == []
    with pytest.raises(ValueError, match="Unknown task guide"):
        workspace.task_guide("conditions")
    assert workspace.publish_lessons() == {"published": 0, "enabled": False}


def test_api_validates_modes_and_serves_free_answers(service):
    app = FastAPI()
    app.include_router(create_scientific_router(service))
    with TestClient(app, base_url="http://127.0.0.1") as client:
        config = client.get("/api/v1/scientific/config").json()
        assert [item["mode"] for item in config["workspace_modes"]] == [
            "pure_agent", "tools_only", "normal", "tools_formatting", "tools_guidance",
        ]
        assert config["default_workspace_mode"] == "normal"
        assert 'id="workspace-mode"' in client.get("/scientific").text
        headers = {"x-scientific-token": config["token"]}
        endpoint = "/api/v1/scientific/turns"
        assert client.post(endpoint, json={"question": "Q", "mode": "bad"}, headers=headers).status_code == 422
        response = client.post(endpoint, json={"question": "Q", "mode": "pure_agent"}, headers=headers)
        assert response.status_code == 202
        identity = response.json()["conversation_id"]
        assert finish(service, identity)["status"] == "completed"
        data = client.get(f"/api/v1/scientific/conversations/{identity}").json()
        assert "Unstructured text" in data["turns"][0]["answer_presentation"]["html"]
        assert client.post(endpoint, json={"question": "Q", "conversation_id": identity, "mode": "normal"},
                           headers=headers).status_code == 422


@pytest.mark.parametrize("mode", list(WorkspaceMode))
@pytest.mark.parametrize("thread", [None, "saved-thread"])
def test_runtime_command_treatments_preserve_model_and_search(tmp_path, monkeypatch, mode, thread):
    runtime = object.__new__(CodexRuntime)
    runtime.executable, runtime.model = "codex", "fixed-model"
    runtime.settings = resolve_research_profile("research")
    monkeypatch.setattr(runtime, "_pure_tool_overrides", lambda _: ["-c", "features.shell_tool=false"])
    command = runtime.command(tmp_path, tmp_path, thread, mode)
    assert ("--output-schema" in command) == (mode in {WorkspaceMode.NORMAL, WorkspaceMode.TOOLS_FORMATTING})
    assert command[command.index("--model") + 1] == "fixed-model"
    assert 'web_search="live"' in command
    assert ('features.shell_tool=false' in command) == (mode == WorkspaceMode.PURE)
    assert ('project_doc_max_bytes=0' in command) == (mode not in {WorkspaceMode.NORMAL, WorkspaceMode.TOOLS_GUIDANCE})
    assert ("resume" in command) == bool(thread)


def test_pure_overrides_disable_each_inherited_mcp_and_fail_closed(tmp_path, monkeypatch):
    runtime = object.__new__(CodexRuntime)
    runtime.executable = "codex"
    features = "\n".join(f"{name} stable true" for name in (
        "shell_tool", "apps", "plugins", "skip_host_skill_discovery", "multi_agent", "browser_use",
        "view_image", "code_mode", "code_mode_host"))
    def command(argv, **kwargs):
        return SimpleNamespace(stdout=features if "features" in argv else json.dumps([
            {"name": "local_data", "enabled": 'mcp_servers.local_data.enabled=false' not in argv},
        ]))
    monkeypatch.setattr("subprocess.run", command)
    overrides = runtime._pure_tool_overrides(tmp_path)
    for item in ('features.shell_tool=false', 'features.apps=false', 'features.plugins=false',
                 'features.multi_agent=false', 'features.view_image=false', 'mcp_servers.local_data.enabled=false',
                 'features.code_mode_host=true'):
        assert item in overrides
    assert 'features.code_mode_host=false' not in overrides
    assert 'features.code_mode=false' not in overrides
    assert 'tools.view_image=false' not in overrides
    for thread in (None, 'saved-thread'):
        runtime.model = None
        runtime.settings = resolve_research_profile('research')
        arguments = runtime.command(tmp_path, tmp_path, thread, WorkspaceMode.PURE)
        assert 'web_search="live"' in arguments
        assert 'features.code_mode_host=true' in arguments
        assert 'features.shell_tool=false' in arguments
    monkeypatch.setattr("subprocess.run", lambda argv, **kwargs: SimpleNamespace(
        stdout=features if "features" in argv else json.dumps([{"name": "local_data", "enabled": True}]),
    ))
    with pytest.raises(RuntimeError, match="could not disable"):
        runtime._pure_tool_overrides(tmp_path)
    features = "shell_tool stable true"
    with pytest.raises(RuntimeError, match="requires a current Codex"):
        runtime._pure_tool_overrides(tmp_path)


@pytest.mark.parametrize("mode", [WorkspaceMode.PURE, WorkspaceMode.TOOLS, WorkspaceMode.TOOLS_GUIDANCE])
def test_empty_free_answer_fails_without_format_repair(service, mode):
    class EmptyRuntime(FreeRuntime):
        def run(self, **kwargs):
            self.calls.append(kwargs)
            return AgentResult({"answer_markdown": ""}, "thread", {})

    service.runtime = EmptyRuntime()
    identity = service.submit("Question", mode=mode)["conversation_id"]
    turn = finish(service, identity)
    assert turn["status"] == "failed"
    assert "answer" not in turn and "repair_attempts" not in turn
    assert len(service.runtime.calls) == 1


@pytest.mark.parametrize("mode", [WorkspaceMode.PURE, WorkspaceMode.TOOLS, WorkspaceMode.TOOLS_GUIDANCE])
def test_real_process_free_text_has_no_schema_or_handoff(tmp_path, monkeypatch, mode):
    attempt = tmp_path / "turns" / "attempt"
    attempt.mkdir(parents=True)
    script = tmp_path / "fake_cli.py"
    script.write_text(
        "import json, os, pathlib, sys\n"
        "attempt = pathlib.Path(sys.argv[1])\n"
        "(attempt / 'agent-final.txt').write_text('Plain final answer.\\n', encoding='utf-8')\n"
        "(attempt / 'environment.json').write_text(json.dumps({'cwd': os.getcwd(), 'path': os.environ.get('PYTHONPATH')}))\n"
        "print(json.dumps({'type': 'thread.started', 'thread_id': 'thread'}))\n"
        "print(json.dumps({'type': 'turn.completed', 'usage': {}}))\n", encoding="utf-8")
    runtime = object.__new__(CodexRuntime)
    runtime.timeout_seconds, runtime.model = 10, None
    runtime.command = lambda *args: [sys.executable, str(script), str(attempt)]
    monkeypatch.setenv("PYTHONPATH", "inherited-repository")
    result = runtime.run(prompt="Just the question", workspace=tmp_path, turn_directory=attempt,
                         thread_id=None, cancel=Event(), on_event=lambda _: None, mode=mode)
    assert result.answer["answer_markdown"] == (attempt / "agent-final.txt").read_bytes().decode("utf-8")
    assert (attempt / "prompt.txt").read_text("utf-8") == "Just the question"
    assert not list(attempt.glob("*schema.json"))
    environment = json.loads((attempt / "environment.json").read_text())
    if mode == WorkspaceMode.PURE:
        assert not Path(environment["cwd"]).is_relative_to(tmp_path)
        assert environment["path"] is None
    else:
        assert Path(environment["cwd"]) == tmp_path
        assert "inherited-repository" in environment["path"]


def test_free_text_rejects_stale_and_linked_output(tmp_path):
    output = tmp_path / "agent-final.txt"
    output.write_text("Old answer", "utf-8")
    with pytest.raises(RuntimeError, match="already contains output"):
        prepare_answer_handoff(tmp_path, tmp_path)
    (tmp_path / "other.txt").hardlink_to(output)
    with pytest.raises(ValueError, match="not a link"):
        load_text_answer(tmp_path, tmp_path)


def test_pure_environment_does_not_mutate_parent(tmp_path):
    original = {"PYTHONPATH": str(tmp_path), "PYTHONSTARTUP": "injected.py", "PATH": "tools"}
    child = runtime_environment(None, search_tool={}, environment=original)
    assert "PYTHONPATH" not in child and "PYTHONSTARTUP" not in child
    assert original["PYTHONPATH"] == str(tmp_path)


class StructuredRuntime(FreeRuntime):
    """Submit a valid conceptual answer without invented scientific evidence."""

    def run(self, **kwargs):
        self.calls.append(kwargs)
        return AgentResult({
            "schema_version": "scientific_answer.v2", "answer_markdown": "A conceptual answer.",
            "sources": [], "molecules": [], "target_molecule_ids": [],
            "steps": [], "routes": [], "claims": [], "evidence_refs": [],
            "uncertainties": [], "needs_user_input": False,
        }, "thread", {"input_tokens": 3})


@pytest.mark.parametrize("mode,guidance,formatting", [
    (WorkspaceMode.TOOLS, False, False),
    (WorkspaceMode.TOOLS_FORMATTING, False, True),
    (WorkspaceMode.TOOLS_GUIDANCE, True, False),
    (WorkspaceMode.NORMAL, True, True),
])
def test_guidance_and_formatting_are_independent_in_api_prompt_and_snapshot(service, mode, guidance, formatting):
    service.runtime = StructuredRuntime() if formatting else FreeRuntime()
    app = FastAPI()
    app.include_router(create_scientific_router(service))
    with TestClient(app, base_url="http://127.0.0.1") as client:
        config = client.get("/api/v1/scientific/config").json()
        html = client.get("/scientific").text
        assert f'value="{mode.value}"' in html
        response = client.post("/api/v1/scientific/turns", json={"question": "Explain", "mode": mode.value},
                               headers={"x-scientific-token": config["token"]})
        assert response.status_code == 202
        identity = response.json()["conversation_id"]
        turn = finish(service, identity)
        assert turn["status"] == "completed", turn.get("error")
        saved = client.get(f"/api/v1/scientific/conversations/{identity}").json()
        assert saved["workspace_mode"] == mode.value
        assert saved["turns"][0]["answer_presentation"]["html"]

    policy = turn["mode_policy"]
    assert policy["task_guidance"] is guidance
    assert policy["structured_answer"] is formatting
    assert turn["answer"]["mode_policy"] == policy == mode_policy(mode)
    prompt = service.runtime.calls[0]["prompt"]
    assert "@DISABLED_RESOURCES@" not in prompt
    assert ("w.task_guide(" in prompt) is guidance
    assert ("w.recall_lessons(" in prompt) is guidance
    assert ("scientific_answer.v2" in prompt) is formatting
    assert ("Default answer presentation" in prompt) is formatting
    assert ("w.finalize_answer(" in prompt) is formatting
    workspace = ScientificWorkspace(service.root / identity)
    context = current_application_context(workspace.store)
    assert context["learning_context"]["enabled"] is guidance
    assert workspace.store.manifest["baseline"]["learning_context"]["enabled"] is guidance
    assert bool(context["learning_context"]["guides"]) is guidance
    resources = context["resources"]
    assert ("agent_instructions/answer_authoring.md" in resources) is formatting
    assert ("presentation/default.md" in resources) is formatting
    assert ("agent_instructions/core.md" in resources) is guidance
    if guidance:
        assert workspace.task_guide("conditions")["text"]
    else:
        assert context["learning_context"] == disabled_learning_context()
        assert workspace.publish_lessons() == {"published": 0, "enabled": False}


@pytest.mark.parametrize("recover", [True, False])
def test_formatting_repair_does_not_enable_task_guidance(service, recover):
    class RepairRuntime(StructuredRuntime):
        def run(self, **kwargs):
            result = super().run(**kwargs)
            if not recover or len(self.calls) == 1:
                result.answer["evidence_refs"] = ["sha256:" + "0" * 64]
            return result

    service.runtime = RepairRuntime()
    identity = service.submit("Explain", mode="tools_formatting")["conversation_id"]
    turn = finish(service, identity)
    assert turn["status"] == ("completed" if recover else "failed"), turn.get("error")
    assert turn["repair_attempts"] == 1
    assert len(service.runtime.calls) == 2
    for call in service.runtime.calls:
        assert call["mode"] == WorkspaceMode.TOOLS_FORMATTING
        assert "w.task_guide(" not in call["prompt"]
        assert "w.recall_lessons(" not in call["prompt"]
        assert "scientific_answer.v2" in call["prompt"]
    assert service.runtime.calls[-1]["thread_id"] == "thread"
    assert "previous submission was rejected" in service.runtime.calls[-1]["prompt"]
    if recover:
        service.submit("Follow up", identity)
        second = finish(service, identity)
        assert second["status"] == "completed", second.get("error")
        assert second["thread_resume"] == "resumed"
        assert second["workspace_mode"] == "tools_formatting"
