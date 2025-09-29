from __future__ import annotations

import io
import json
from typing import Any, Dict, Iterable, List

import pytest

from clients.mcp_rules_client import McpRulesClient


class FakeStdout:
    def __init__(self, responses: Iterable[Dict[str, Any]]) -> None:
        self._responses = [json.dumps(resp) + "\n" for resp in responses]

    def readline(self) -> str:
        if not self._responses:
            return ""
        return self._responses.pop(0)


class FakeStdin:
    def __init__(self) -> None:
        self.buffer: List[str] = []

    def write(self, data: str) -> None:
        self.buffer.append(data)

    def flush(self) -> None:  # pragma: no cover - noop
        return None


class FakeProcess:
    def __init__(self, responses: Iterable[Dict[str, Any]]) -> None:
        self.stdin = FakeStdin()
        self.stdout = FakeStdout(responses)
        self.stderr = io.StringIO()
        self._alive = True

    def poll(self) -> int | None:
        return None if self._alive else 0

    def terminate(self) -> None:  # pragma: no cover - noop
        self._alive = False

    def wait(self, timeout: float | None = None) -> None:  # pragma: no cover
        self._alive = False

    def kill(self) -> None:  # pragma: no cover - noop
        self._alive = False


@pytest.fixture
def patched_popen(monkeypatch: pytest.MonkeyPatch) -> None:
    responses = [
        {"jsonrpc": "2.0", "id": 1, "result": {"status": "ok"}},
        {
            "jsonrpc": "2.0",
            "id": 2,
            "result": {
                "matched_playbooks": ["PB-1"],
                "suggestions": [{"playbook_id": "PB-1"}],
            },
        },
    ]

    def fake_popen(*_args: Any, **_kwargs: Any) -> FakeProcess:
        return FakeProcess(responses)

    monkeypatch.setattr("clients.mcp_rules_client.subprocess.Popen", fake_popen)


def test_client_formats_requests(monkeypatch: pytest.MonkeyPatch, patched_popen: None) -> None:
    client = McpRulesClient("python fake_server.py", "dummy.json", startup_timeout_s=0.1)
    client.start()
    result = client.apply("A>>B", {"electrophile": {"class": "aryl chloride"}}, max_suggestions=1)
    assert result["matched_playbooks"] == ["PB-1"]


def test_client_detects_dead_process(monkeypatch: pytest.MonkeyPatch) -> None:
    class DeadProcess(FakeProcess):
        def __init__(self) -> None:
            super().__init__(responses=[])
            self._alive = False

    monkeypatch.setattr("clients.mcp_rules_client.subprocess.Popen", lambda *args, **kwargs: DeadProcess())
    client = McpRulesClient("python fake_server.py", "dummy.json", startup_timeout_s=0.05)
    with pytest.raises(TimeoutError):
        client.start()