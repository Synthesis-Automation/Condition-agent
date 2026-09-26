"""Conversation harness, persisted evidence, runtime failure, and local API regressions."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from threading import Event
from time import monotonic, sleep
from typing import Any

from fastapi.testclient import TestClient
import pytest

from app.web_api.main import create_app
from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.agent_runtime import AgentResult, AgentStopped, CodexRuntime
from chem_coworker.scientific_workspace.baseline import code_manifest, environment_versions
from chem_coworker.scientific_workspace.conversation import ConversationService


ROOT = Path(__file__).resolve().parents[2]


class RecordedRuntime:
    """Test double for the model only; execute real recorded reaction chemistry."""

    def __init__(self, *, missing_reference: bool = False, wait: bool = False) -> None:
        self.threads: list[str | None] = []
        self.missing_reference = missing_reference
        self.wait = wait

    def describe(self) -> dict[str, Any]:
        return {"runtime": "test_double", "version": "1", "model": "none"}

    def run(self, *, prompt, workspace, turn_directory, thread_id, cancel, on_event):
        self.threads.append(thread_id)
        assert "Do not invent" in prompt
        on_event({"type": "thread.started", "thread_id": "test-thread"})
        if self.wait:
            assert cancel.wait(10), "Test failed to cancel the active turn"
            raise AgentStopped("cancelled")
        scientific = ScientificWorkspace(workspace)
        event = scientific.run("analyze_reaction", {"reaction_smiles": "CCBr.N>>CCN"})
        reference = "sha256:" + "0" * 64 if self.missing_reference else event.artifact_ref
        return AgentResult({
            "schema_version": "scientific_answer.v2",
            "sources": [], "molecules": [], "target_molecule_ids": [],
            "steps": [], "routes": [], "claims": [],
            "answer_markdown": "Recorded graph interpretation; this is not experimental validation.",
            "evidence_refs": [reference], "uncertainties": ["No experimental outcome verified"],
            "needs_user_input": False,
        }, "test-thread", {"input_tokens": 10})


@pytest.fixture
def service(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    # Avoid rehashing GB datasets or repeating registry audit in harness-only tests.
    monkeypatch.setattr("chem_coworker.scientific_workspace.baseline.capture_baseline", lambda *_: {
        "repository": str(ROOT), "code_files": code_manifest(ROOT),
        "environment": environment_versions(), "artifacts": {},
        "validation_status": "development_snapshot_not_release_validated",
    })
    instance = ConversationService(tmp_path / "conversations", runtime=RecordedRuntime())
    yield instance
    instance.close()


def finish(service: ConversationService, identity: str) -> dict[str, Any]:
    deadline = monotonic() + 20
    while monotonic() < deadline:
        result = service.get(identity)["turns"][-1]
        if result["status"] not in {"queued", "preparing", "running"} and service._active is None:
            return result
        sleep(0.02)
    raise AssertionError("Conversation did not finish")


def test_real_evidence_answer_and_follow_up_survive_service_restart(service: ConversationService) -> None:
    first = service.submit("Analyze CCBr.N>>CCN")
    turn = finish(service, first["conversation_id"])
    assert turn["status"] == "completed", turn
    assert turn["answer"]["review_status"] == "unreviewed"
    evidence = service.artifact(first["conversation_id"], turn["answer"]["evidence_refs"][0])
    assert evidence["operation"] == "analyze_reaction"
    assert evidence["execution_status"] == "completed"
    reopened = ConversationService(service.root, runtime=service.runtime)
    try:
        reopened.submit("What remains unknown?", first["conversation_id"])
        assert finish(reopened, first["conversation_id"])["status"] == "completed"
        assert service.runtime.threads == [None, "test-thread"]
        assert len(reopened.get(first["conversation_id"])["turns"]) == 2
    finally:
        reopened.close()


def test_invented_citation_is_a_failed_turn_not_a_supported_answer(service: ConversationService) -> None:
    service.runtime = RecordedRuntime(missing_reference=True)
    identity = service.submit("Analyze a reaction")["conversation_id"]
    turn = finish(service, identity)
    assert turn["status"] == "failed"
    assert "answer" not in turn
    assert len(service.runtime.threads) == 2
    assert any(event.kind == "call" for event in ScientificWorkspace(service.root / identity).store.events())


def test_explicit_runtime_setting_change_starts_new_thread_with_saved_history(service: ConversationService) -> None:
    identity = service.submit("Analyze CCBr.N>>CCN")["conversation_id"]
    first = finish(service, identity)
    assert first["runtime_requested"]["model"] == "none"
    assert first["answer"]["self_review_status"] == "not_recorded_for_final_draft"
    runtime = RecordedRuntime()
    runtime.describe = lambda: {"runtime": "test_double", "version": "1", "model": "changed"}
    service.runtime = runtime
    service.submit("Continue with the saved scientific evidence", identity)
    second = finish(service, identity)
    assert second["status"] == "completed", second
    assert second["thread_resume"] == "configuration_changed_new_thread"
    assert runtime.threads == [None]
    assert len(ScientificWorkspace(service.root / identity).store.events()) >= 6
    assert (service.root / identity / "turns" / second["id"] / "capabilities.json").is_file()


def test_call_summary_preserves_result_and_measures_execution_phases(service: ConversationService) -> None:
    identity = service.submit("Analyze CCBr.N>>CCN")["conversation_id"]
    finish(service, identity)
    workspace = ScientificWorkspace(service.root / identity)
    event = next(event for event in workspace.store.events() if event.kind == "call")
    payload = workspace.store.read_artifact(event.artifact_ref)
    summary = workspace.call_summary(event)
    assert payload["result_bytes"] > 0
    assert set(payload["timings"]) == {"baseline_and_evidence_seconds", "operation_seconds", "serialization_seconds"}
    assert all(value >= 0 for value in payload["timings"].values())
    assert summary["event"]["artifact_ref"] == event.artifact_ref
    assert workspace.store.read_artifact(event.artifact_ref) == payload


def test_one_repair_uses_same_thread_and_preserves_rejected_answer(service: ConversationService) -> None:
    class RepairRuntime(RecordedRuntime):
        def run(self, **kwargs):
            self.missing_reference = not self.threads
            return super().run(**kwargs)

    service.runtime = RepairRuntime()
    identity = service.submit("Analyze a reaction")['conversation_id']
    turn = finish(service, identity)
    assert turn["status"] == "completed", turn
    assert turn["repair_attempts"] == 1
    assert service.runtime.threads == [None, "test-thread"]
    assert len(turn["answer"]["attempt_usage"]) == 2
    assert (service.root / identity / "turns" / turn["id"] / "rejected-answer.json").is_file()
    events = ScientificWorkspace(service.root / identity).store.events()
    assert sum(event.kind == "agent_answer_rejected" for event in events) == 1


def test_cancellation_preserves_conversation_and_prevents_concurrent_turns(service: ConversationService) -> None:
    service.runtime = RecordedRuntime(wait=True)
    identity = service.submit("Wait for additional evidence")["conversation_id"]
    with pytest.raises(RuntimeError, match="running"):
        service.submit("Second question")
    assert service.cancel(identity)
    assert finish(service, identity)["status"] == "cancelled"
    assert not (service.root / identity / ".conversation.lock").exists()


def test_baseline_change_blocks_follow_up_before_model_execution(service: ConversationService, monkeypatch) -> None:
    identity = service.submit("Analyze CCBr.N>>CCN")["conversation_id"]
    finish(service, identity)
    monkeypatch.setattr("chem_coworker.scientific_workspace.baseline.code_manifest", lambda _: {})
    service.submit("Continue", identity)
    turn = finish(service, identity)
    assert turn["status"] == "failed"
    assert "changed" in turn["error"]["message"]
    assert service.runtime.threads == [None]


def test_ids_questions_and_cross_process_lock_are_checked(service: ConversationService) -> None:
    for question in (" ", "a" * 20001):
        with pytest.raises(ValueError):
            service.submit(question)
    with pytest.raises(ValueError):
        service.get("../../outside")
    with pytest.raises(FileNotFoundError):
        service.submit("Question", "0" * 32)
    identity = service.submit("Question")["conversation_id"]
    finish(service, identity)
    lock = service.root / identity / ".conversation.lock"
    lock.write_text("unknown worker", "utf-8")
    with pytest.raises(RuntimeError, match="locked"):
        service.submit("Continue", identity)


def test_interrupted_worker_is_visible_without_removing_its_lock(service: ConversationService) -> None:
    identity = service.submit("Question")["conversation_id"]
    turn = finish(service, identity)
    path = service.root / identity / "turns" / turn["id"] / "turn.json"
    turn["status"] = "running"
    path.write_text(json.dumps(turn), "utf-8")
    lock = service.root / identity / ".conversation.lock"
    lock.write_text("123", "utf-8")
    assert service.get(identity)["turns"][-1]["status"] == "interrupted"
    assert lock.exists()


def test_atomic_json_retries_reader_sharing_violation(tmp_path: Path, monkeypatch) -> None:
    from chem_coworker.scientific_workspace import store

    original = store.os.replace
    attempts = []

    def replace(source, destination):
        attempts.append(source)
        if len(attempts) == 1:
            raise PermissionError("Windows reader currently has the destination open")
        return original(source, destination)

    monkeypatch.setattr(store.os, "replace", replace)
    path = tmp_path / "state.json"
    store._write_json(path, {"status": "completed"})
    assert json.loads(path.read_text("utf-8")) == {"status": "completed"}
    assert len(attempts) == 2
    assert not list(tmp_path.glob(".writing-*"))


def test_json_reader_retries_sharing_violation_but_not_invalid_json(tmp_path: Path, monkeypatch) -> None:
    from chem_coworker.scientific_workspace.store import _read_json

    path = tmp_path / "state.json"
    path.write_text('{"status":"completed"}', "utf-8")
    original = Path.read_text
    attempts = []

    def read(instance, *args, **kwargs):
        attempts.append(instance)
        if len(attempts) == 1:
            raise PermissionError("Transient replacement sharing violation")
        return original(instance, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", read)
    assert _read_json(path) == {"status": "completed"}
    assert len(attempts) == 2
    path.write_text("broken", "utf-8")
    with pytest.raises(json.JSONDecodeError):
        _read_json(path)
    assert len(attempts) == 3


def test_optional_web_profile_origin_token_and_evidence(service: ConversationService) -> None:
    client = TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False), base_url="http://127.0.0.1")
    assert client.get("/scientific").status_code == 200
    assert "Ask a chemistry question" in client.get("/scientific").text
    token = client.get("/api/v1/scientific/config").json()["token"]
    endpoint = "/api/v1/scientific/turns"
    assert client.post(endpoint, json={"question": "Hello"}).status_code == 403
    headers = {"x-scientific-token": token}
    assert client.post(endpoint, json={"question": "Hello"}, headers={**headers, "origin": "https://attacker.example"}).status_code == 403
    assert client.get("/api/v1/scientific/config", headers={"host": "attacker.example"}).status_code == 403
    response = client.post(endpoint, json={"question": "Analyze CCBr.N>>CCN"}, headers=headers)
    assert response.status_code == 202
    identity = response.json()["conversation_id"]
    turn = finish(service, identity)
    ref = turn["answer"]["evidence_refs"][0]
    assert client.get(f"/api/v1/scientific/conversations/{identity}/artifacts/{ref}").json()["operation"] == "analyze_reaction"
    assert client.get("/api/v1/scientific/conversations/invalid").status_code == 422
    assert client.post(endpoint, json={"question": "Q", "executable": "arbitrary"}, headers=headers).status_code == 422
    focused = TestClient(create_app(runtime=object(), scientific_service=service), base_url="http://127.0.0.1")
    assert focused.get("/scientific").status_code == 404
    assert focused.get("/api/v1/scientific/config").status_code == 404


def test_web_activity_finds_and_cancels_the_owned_investigation(service: ConversationService) -> None:
    service.runtime = RecordedRuntime(wait=True)
    client = TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False),
                        base_url="http://127.0.0.1")
    endpoint = "/api/v1/scientific/activity"
    assert client.get(endpoint).json() == {"active": None}
    token = client.get("/api/v1/scientific/config").json()["token"]
    submitted = service.submit("Investigate this reaction")
    active = client.get(endpoint).json()["active"]
    assert active["conversation_id"] == submitted["conversation_id"]
    assert active["id"] == submitted["turn_id"]
    assert active["status"] in {"queued", "preparing", "running"}
    assert active["progress"] == []
    assert "answer" not in active
    assert client.get(endpoint, headers={"origin": "https://attacker.example"}).status_code == 403
    response = client.post(f"/api/v1/scientific/conversations/{active['conversation_id']}/cancel",
                           headers={"x-scientific-token": token})
    assert response.json() == {"cancellation_requested": True}
    assert finish(service, active["conversation_id"])["status"] == "cancelled"
    assert client.get(endpoint).json() == {"active": None}


def test_web_activity_does_not_adopt_an_interrupted_turn(service: ConversationService) -> None:
    identity = service.submit("Question")["conversation_id"]
    turn = finish(service, identity)
    path = service.root / identity / "turns" / turn["id"] / "turn.json"
    turn["status"] = "running"
    path.write_text(json.dumps(turn), "utf-8")
    client = TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False),
                        base_url="http://127.0.0.1")
    assert client.get("/api/v1/scientific/activity").json() == {"active": None}


def test_runtime_command_never_uses_shell_or_last_global_session(tmp_path: Path) -> None:
    runtime = object.__new__(CodexRuntime)
    runtime.executable = "codex"
    runtime.model = None
    command = runtime.command(tmp_path, tmp_path / "turn", "specific-thread")
    assert command[:4] == ["codex", "exec", "--sandbox", "workspace-write"]
    assert command[command.index("resume") + 1] == "specific-thread"
    assert "--last" not in command
    assert "--dangerously-bypass-approvals-and-sandbox" not in command
    assert command[-1] == "-"


@pytest.mark.parametrize("mode", ["success", "failed", "no_final", "timeout"])
def test_native_process_events_completion_and_timeout(tmp_path: Path, mode: str) -> None:
    script = tmp_path / "fake_runtime.py"
    script.write_text('''import json, pathlib, sys, time
mode, output = sys.argv[1:]
print(json.dumps({"type":"thread.started","thread_id":"native-thread"}), flush=True)
if mode == "timeout":
    time.sleep(30)
elif mode == "failed":
    print(json.dumps({"type":"turn.failed","error":{"message":"provider unavailable"}}), flush=True)
else:
    if mode != "no_final":
        pathlib.Path(output).write_text(json.dumps({"answer_markdown":"Test", "evidence_refs":[], "uncertainties":[], "needs_user_input":False}), encoding="utf-8")
    print(json.dumps({"type":"turn.completed","usage":{"output_tokens":4}}), flush=True)
''', "utf-8")
    runtime = object.__new__(CodexRuntime)
    runtime.timeout_seconds = 1 if mode == "timeout" else 10
    runtime.command = lambda *_: [sys.executable, str(script), mode, str(tmp_path / "agent-final.json")]
    events = []
    kwargs = dict(prompt="用户问题; $(not shell code)", workspace=tmp_path, turn_directory=tmp_path,
                  thread_id=None, cancel=Event(), on_event=events.append)
    if mode == "success":
        result = runtime.run(**kwargs)
        assert result.thread_id == "native-thread"
        assert result.usage == {"output_tokens": 4}
        assert len(events) == 2
        assert (tmp_path / "prompt.txt").read_text("utf-8") == kwargs["prompt"]
    elif mode == "timeout":
        with pytest.raises(AgentStopped, match="timed_out"):
            runtime.run(**kwargs)
    else:
        with pytest.raises(RuntimeError):
            runtime.run(**kwargs)
