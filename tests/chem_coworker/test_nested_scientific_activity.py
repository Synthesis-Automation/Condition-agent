"""Recorded scientific failures remain visible independently of outer shell success."""

from __future__ import annotations

import json
from pathlib import Path
from threading import Event
from time import monotonic, sleep

import pytest

from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.runtime.activity import (
    ActivityHistory,
    ScientificActivityCursor,
)
from chem_coworker.scientific_workspace.runtime.agent_runtime import AgentResult
from chem_coworker.scientific_workspace.core.baseline import environment_versions
from chem_coworker.scientific_workspace.runtime.conversation import ConversationService
from chem_coworker.scientific_workspace.core.store import InvestigationStore


def test_cursor_reads_only_current_turn_and_deduplicates_events(tmp_path: Path) -> None:
    store = InvestigationStore.create(tmp_path / "run", objective="Check failures", baseline={})
    old = store.append("call", {"operation": "old_call", "execution_status": "error"})
    start = store.append("user_message", {"text": "New turn"})
    history = ActivityHistory()
    cursor = ScientificActivityCursor(store, after_sequence=start.sequence)
    failure = store.append("call", {
        "operation": "assess_route_proposal", "execution_status": "error",
        "error": {"type": "ValueError", "message": "Invalid evidence reference"},
        "arguments": {"source_ref": "missing"}, "duration_seconds": 0.1,
    })
    success = store.append("call", {
        "operation": "disconnect_target", "execution_status": "completed",
        "arguments": {"target_smiles": "CCN"},
        "result": {"strategies": [{"id": "one"}], "valid": False},
    })
    rows = cursor.drain(history)
    assert [row["event_sequence"] for row in rows] == [failure.sequence, success.sequence]
    assert rows[0]["status"] == "failed"
    assert rows[0]["failure_detail"] == "ValueError: Invalid evidence reference"
    assert rows[1]["status"] == "completed"  # Scientific uncertainty is not execution failure.
    assert "target_smiles: CCN" in rows[1]["detail"]
    assert "1 strategies returned" in rows[1]["detail"]
    assert old.artifact_ref not in json.dumps(rows)
    assert cursor.drain(history) == []
    assert len(history.rows) == 2


def test_literature_failure_capture_and_extraction_are_distinguished(tmp_path: Path) -> None:
    store = InvestigationStore.create(tmp_path / "run", objective="Read sources", baseline={})
    cursor = ScientificActivityCursor(store, after_sequence=0)
    history = ActivityHistory()
    store.append("literature_source", {
        "source_url": "https://example.org/paper", "retrieval_status": "failed",
        "extraction": {"status": "not_attempted"},
        "error": {"type": "PermissionError", "message": "Socket denied"},
    })
    store.append("literature_source", {
        "source_url": "https://example.org/paper", "retrieval_status": "not_performed",
        "extraction": {"status": "agent_supplied", "text": "Long source text excluded"},
    })
    store.append("literature_source", {
        "source_url": "https://example.org/paper.pdf", "retrieval_status": "completed",
        "extraction": {"status": "failed", "error": "Invalid PDF document"},
    })
    rows = cursor.drain(history)
    assert [row["status"] for row in rows] == ["failed", "completed", "failed"]
    assert rows[0]["failure_detail"] == "PermissionError: Socket denied"
    assert "not HTTP-verified" in rows[1]["detail"]
    assert "https://example.org/paper" in rows[1]["detail"]
    assert rows[2]["title"] == "Extract literature text"
    assert rows[2]["failure_detail"] == "Invalid PDF document"
    assert "Long source text excluded" not in json.dumps(rows)


@pytest.mark.parametrize("terminal", ["completed", "failed"])
def test_turn_logs_nested_errors_on_heartbeat_and_finalization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, terminal: str,
) -> None:
    ready, release = Event(), Event()
    monkeypatch.setattr("chem_coworker.scientific_workspace.core.baseline.capture_baseline", lambda *_, **kwargs: {
        "repository": str(tmp_path), "artifacts": {}, "code_files": {},
        "environment": environment_versions(),
    })
    monkeypatch.setattr("chem_coworker.scientific_workspace.runtime.conversation.verify_baseline", lambda *_: None)

    class Runtime:
        def describe(self):
            return {"runtime": "nested_activity_test"}

        def run(self, *, workspace, on_event, **kwargs):
            store = ScientificWorkspace(workspace).store
            store.append("call", {
                "operation": "assess_route_proposal", "execution_status": "error",
                "error": {"type": "ValueError", "message": "Invalid evidence reference"},
            })
            on_event({"type": "runtime.heartbeat"})
            on_event({"type": "runtime.heartbeat"})
            ready.set()
            assert release.wait(10)
            # No runtime events follow these writes: finalization must collect them.
            store.append("literature_source", {
                "url": "https://example.org/paper", "retrieval_status": "failed",
                "error": {"type": "PermissionError", "message": "Socket denied"},
            })
            if terminal == "failed":
                raise RuntimeError("Provider disconnected")
            return AgentResult({
                "schema_version": "scientific_answer.v2", "sources": [], "molecules": [],
                "target_molecule_ids": [], "steps": [], "routes": [], "claims": [],
                "answer_markdown": "The attempted checks failed; no scientific claim was verified.",
                "evidence_refs": [], "uncertainties": ["Evidence missing"], "needs_user_input": False,
            }, "nested-test", {})

    service = ConversationService(tmp_path / "conversations", runtime=Runtime(), repository=tmp_path)
    try:
        submitted = service.submit("Check a route")
        identity, turn_id = submitted["conversation_id"], submitted["turn_id"]
        assert ready.wait(10)
        live = service.get(identity)["turns"][0]
        assert live["status"] == "running"
        assert len(live["progress"]) == 1
        assert live["progress"][0]["failure_detail"] == "ValueError: Invalid evidence reference"
        release.set()
        deadline = monotonic() + 10
        while service._active is not None and monotonic() < deadline:
            sleep(0.02)
        final = service.get(identity)["turns"][0]
        assert final["status"] == terminal, final
        assert len(final["progress"]) == 2
        records = [json.loads(line) for line in service.debug_log(identity, turn_id).splitlines()]
        recorded = [record for record in records if "event_sequence" in record]
        assert len(recorded) == 2
        assert all(record["artifact_ref"] == record["activity"]["artifact_ref"] for record in recorded)
        assert all(record["activity"]["failure_detail"] for record in recorded)
        assert not any(record["kind"] == "agent_update" for record in records)
    finally:
        release.set()
        service.close()
