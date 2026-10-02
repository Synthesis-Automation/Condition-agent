"""Local tool discovery, child environment, and blocked literature recovery."""

from __future__ import annotations

import errno
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
from threading import Event

import pytest

from chem_coworker.scientific_workspace.adapters import literature
from chem_coworker.scientific_workspace.runtime import runtime_environment as environment_tools
from chem_coworker.scientific_workspace.runtime.agent_runtime import CodexRuntime
from chem_coworker.scientific_workspace.core.store import InvestigationStore


def test_discovers_and_validates_editor_bundled_ripgrep_when_path_is_missing(tmp_path, monkeypatch):
    bundled = tmp_path / ".vscode/extensions/openai.chatgpt-test/bin/windows-x86_64/rg.exe"
    bundled.parent.mkdir(parents=True)
    bundled.write_bytes(b"test executable")
    monkeypatch.setattr(environment_tools.shutil, "which", lambda *args, **kwargs: None)
    observed = []

    def run(command, **kwargs):
        observed.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, "ripgrep 14.1.0\nfeatures: test\n", "")

    monkeypatch.setattr(environment_tools.subprocess, "run", run)
    result = environment_tools.discover_ripgrep(environment={"PATH": ""}, home=tmp_path)
    assert result["status"] == "available"
    assert result["path"] == str(bundled.resolve())
    assert result["source"] == "editor_extension"
    assert result["version"] == "ripgrep 14.1.0"
    assert observed[0][0] == [str(bundled.resolve()), "--version"]
    assert "shell" not in observed[0][1]


def test_discovers_native_runtime_sibling_and_rejects_broken_path_binary(tmp_path, monkeypatch):
    broken = tmp_path / "broken/rg.exe"
    broken.parent.mkdir()
    broken.write_bytes(b"broken")
    runtime = tmp_path / "bundle/codex.exe"
    runtime.parent.mkdir()
    runtime.write_bytes(b"test runtime")
    bundled = runtime.parent / ("rg.exe" if os.name == "nt" else "rg")
    bundled.write_bytes(b"test search")
    monkeypatch.setattr(environment_tools.shutil, "which", lambda *args, **kwargs: str(broken))
    checked = []

    def run(command, **kwargs):
        checked.append(command[0])
        return subprocess.CompletedProcess(command, 1 if command[0] == str(broken) else 0,
                                           "ripgrep 14.1.0", "")

    monkeypatch.setattr(environment_tools.subprocess, "run", run)
    result = environment_tools.discover_ripgrep(
        environment={"PATH": ""}, home=tmp_path, runtime_executable=str(runtime),
    )
    assert checked == [str(broken.resolve()), str(bundled.resolve())]
    assert result["source"] == "runtime_bundle"


def test_missing_or_unrecognized_tool_reports_fallback_without_installing(tmp_path, monkeypatch):
    impostor = tmp_path / "rg.exe"
    impostor.write_bytes(b"test")
    monkeypatch.setattr(environment_tools.shutil, "which", lambda *args, **kwargs: str(impostor))
    monkeypatch.setattr(environment_tools.subprocess, "run", lambda command, **kwargs:
                        subprocess.CompletedProcess(command, 0, "something else", ""))
    result = environment_tools.discover_ripgrep(environment={"PATH": ""}, home=tmp_path)
    assert result["status"] == "unavailable"
    assert result["path"] is None
    assert "Select-String" in result["fallback"]


def test_search_path_is_configured_only_for_child_and_preserves_existing_entries(tmp_path):
    directory = tmp_path / "bundled tools"
    existing = tmp_path / "existing"
    parent = {"Path": os.pathsep.join([str(existing), str(directory)]), "PYTHONPATH": "custom"}
    before = dict(parent)
    child = environment_tools.runtime_environment(
        tmp_path / "repo", search_tool={"status": "available", "path": str(directory / "rg.exe")},
        environment=parent,
    )
    assert parent == before
    assert child["Path"].split(os.pathsep) == [str(directory), str(existing)]
    assert "PATH" not in child
    assert child["PYTHONPATH"] == str(tmp_path / "repo") + os.pathsep + "custom"
    assert child["PYTHONDONTWRITEBYTECODE"] == "1"


def test_runtime_child_receives_search_path_and_emits_idle_poll_signal(tmp_path):
    script = tmp_path / "runtime_double.py"
    script.write_text(
        "import json, os, pathlib, sys, time\n"
        "prompt = sys.stdin.read()\n"
        "print(json.dumps({'type':'thread.started', 'thread_id':'test-thread'}), flush=True)\n"
        "time.sleep(1.2)\n"
        "pathlib.Path('answer-draft.json').write_text(json.dumps({'path':os.environ['PATH'], 'prompt':prompt}))\n"
        "pathlib.Path('agent-final.json').write_text(json.dumps({'schema_version':'scientific_answer_handoff.v1', 'answer_file':'answer-draft.json'}))\n"
        "print(json.dumps({'type':'turn.completed'}), flush=True)\n", "utf-8",
    )
    runtime = object.__new__(CodexRuntime)
    runtime.timeout_seconds = 10
    runtime.search_tool = {"status": "available", "path": str(tmp_path / "tools/rg.exe")}
    runtime.command = lambda *_: [sys.executable, str(script)]
    events = []
    parent_path = os.environ.get("PATH", "")
    result = runtime.run(prompt="Original prompt", workspace=tmp_path, turn_directory=tmp_path,
                         thread_id=None, cancel=Event(), on_event=events.append)
    assert result.answer["path"].split(os.pathsep)[0] == str(tmp_path / "tools")
    assert result.answer["prompt"].startswith("Original prompt\n\nANSWER FILE HANDOFF")
    assert "answer-draft.json" in result.answer["prompt"]
    assert os.environ.get("PATH", "") == parent_path
    assert {"type": "runtime.heartbeat"} in events
    request = json.loads((tmp_path / "runtime-request.json").read_text("utf-8"))
    assert request["configuration"]["local_tools"]["rg"]["path"] == runtime.search_tool["path"]


@pytest.mark.parametrize("windows", [False, True])
def test_permission_failure_records_browser_capture_recovery_and_preserves_error(tmp_path, monkeypatch, windows):
    store = InvestigationStore.create(tmp_path / "study", objective="Test recovery", baseline={})
    denied = PermissionError(errno.EACCES, "Socket access denied")
    if windows:
        denied = OSError("An attempt was made to access a socket in a forbidden way")
        denied.winerror = 10013

    def fail(*args, **kwargs):
        raise denied

    monkeypatch.setattr(literature, "_request", fail)
    event = literature.fetch_source(store, "https://example.org/patent")
    source = store.read_artifact(event.artifact_ref)
    assert source["retrieval_status"] == "failed"
    assert source["error"]["category"] == "network_permission_denied"
    assert source["error"]["type"] == type(denied).__name__
    assert source["error"]["message"] == str(denied)
    assert source["snapshot"] is None
    assert source["recovery"]["action"] == "capture_browser_source"
    assert "w.capture_source" in source["recovery"]["instructions"]
    assert "Do not repeat" in source["recovery"]["retry_guidance"]
    inspected = literature.inspect_source(store, event.artifact_ref)
    assert inspected["error"] == source["error"]
    assert inspected["recovery"] == source["recovery"]
    captured = literature.capture_source(store, "Exact browser-obtained passage", url="https://example.org/patent")
    assert store.read_artifact(captured.artifact_ref)["retrieval_status"] == "not_performed"
    assert store.read_artifact(event.artifact_ref) == source


def test_socket_permission_denial_stops_equivalent_address_retries(monkeypatch):
    attempts = []
    closed = []

    class DeniedSocket:
        def settimeout(self, value):
            pass

        def connect(self, address):
            attempts.append(address)
            raise PermissionError(errno.EACCES, "Denied")

        def close(self):
            closed.append(True)

    monkeypatch.setattr(socket, "getaddrinfo", lambda *args, **kwargs: [
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("93.184.216.34", 80)),
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("93.184.216.35", 80)),
    ])
    monkeypatch.setattr(socket, "socket", lambda *args: DeniedSocket())
    with pytest.raises(PermissionError):
        literature._request("http://example.org/source", timeout=1)
    assert attempts == [("93.184.216.34", 80)]
    assert len(closed) == 1
