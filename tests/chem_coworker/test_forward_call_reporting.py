"""An optional forward deadline must remain visible in saved calls, CLI and UI."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.__main__ import main
from chem_coworker.scientific_workspace.runtime.activity import ActivityHistory
from chem_coworker.scientific_workspace.adapters.operations import ScientificOperations


@pytest.mark.parametrize("status", ["timed_out", "error"])
def test_forward_incomplete_result_keeps_diagnostics_and_is_not_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str], status: str,
) -> None:
    store = InvestigationStore.create(tmp_path / "run", objective="Inspect one route step", baseline={})
    source = store.append("call", {"operation": "assess_route_proposal", "execution_status": "completed"})
    result = {
        "schema_version": "route_step_forward_investigation.v1", "execution_status": status,
        "error": {"type": "TimeoutError" if status == "timed_out" else "ValueError",
                  "message": "Optional check did not finish"},
        "assessment": None, "source_ref": source.artifact_ref,
        "experimental_feasibility": "not_established",
        "stages": [{"stage": "load_forward_library", "status": "started"}],
    }
    monkeypatch.setattr("chem_coworker.scientific_workspace.workspace.verify_baseline", lambda *_: None)
    monkeypatch.setattr(ScientificOperations, "invoke", lambda *_: result)
    request = {"source_ref": source.artifact_ref, "step_id": "step-1",
               "question": "Could a competing product change the route?"}
    arguments = tmp_path / "request.json"
    arguments.write_text(json.dumps(request), encoding="utf-8")
    assert main(["run", str(store.root), "assess_route_step_forward", "--input", str(arguments)]) == 1
    assert json.loads(capsys.readouterr().out)["execution_status"] == status
    workspace = ScientificWorkspace(store.root)
    event = list(workspace.store.events())[-1]
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == status
    assert payload["result"] == result
    assert payload["error"] == result["error"]
    assert event.evidence_refs == (source.artifact_ref,)
    row = ActivityHistory().observe_scientific({
        "sequence": event.sequence, "kind": "call", "artifact_ref": event.artifact_ref,
        "created_at": event.created_at,
    }, payload)
    assert row["status"] == "failed"
    assert "Optional check did not finish" in row["failure_detail"]
    assert "step_id: step-1" in row["detail"]
    assert workspace.store.summary()["status"] == "active"


@pytest.mark.parametrize("chemistry_changed", [False, True])
def test_forward_replay_compares_chemistry_despite_new_log_paths_and_timings(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, chemistry_changed: bool,
) -> None:
    store = InvestigationStore.create(tmp_path / "run", objective="Replay a selected step", baseline={})
    result = {
        "execution_status": "completed", "assessment": {"disposition": "supported"},
        "provenance": {"forward_library_sha256": "unchanged"},
        "execution": {
            "timings": {"elapsed_seconds": 1.0, "timeout_seconds": 30},
            "stages": [{"stage": "assess_selected_step", "status": "completed",
                        "elapsed_seconds": 1.0, "duration_seconds": 0.5}],
            "diagnostics": {"stdout.log": "first/stdout.log"},
        },
    }
    source = store.append("call", {"operation": "assess_route_step_forward", "arguments": {},
                                  "execution_status": "completed", "result": result})
    actual = deepcopy(result)
    actual["execution"]["timings"]["elapsed_seconds"] = 2.0
    actual["execution"]["stages"][0].update(elapsed_seconds=2.0, duration_seconds=1.5)
    actual["execution"]["diagnostics"]["stdout.log"] = "second/stdout.log"
    if chemistry_changed:
        actual["assessment"]["disposition"] = "competitive"
    monkeypatch.setattr("chem_coworker.scientific_workspace.workspace.verify_baseline", lambda *_, **__: None)
    monkeypatch.setattr(ScientificOperations, "invoke", lambda *_: actual)
    workspace = ScientificWorkspace(store.root)
    event = workspace.replay(source.artifact_ref)
    replay = workspace.store.read_artifact(event.artifact_ref)
    assert replay["matches"] is not chemistry_changed
    assert replay["result"] == actual
    assert workspace.store.read_artifact(source.artifact_ref)["result"] == result
