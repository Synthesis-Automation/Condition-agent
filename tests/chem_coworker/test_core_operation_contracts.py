"""New capabilities use the same recording and replay core as existing tools."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping

import pytest

from chem_coworker.scientific_workspace import (
    InvestigationStore, OperationDefinition, ScientificWorkspace,
)
from chem_coworker.scientific_workspace.operations import ScientificOperations


class TrialOperations(ScientificOperations):
    """A future capability declared outside the workspace orchestrator."""

    DEFINITIONS = ScientificOperations.DEFINITIONS + (
        OperationDefinition(
            "inspect_trial", contract_version="2", evidence_arguments=("observations",),
            execution_status_field="execution_status",
        ),
    )

    def inspect_trial(
        self, observations: list[str], status: str = "unknown",
        execution_status: str = "completed",
    ) -> dict[str, Any]:
        """Keep scientific uncertainty independent of execution completion."""
        return {"observations": observations, "status": status,
                "execution_status": execution_status}


@pytest.fixture
def workspace(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> ScientificWorkspace:
    store = InvestigationStore.create(
        tmp_path / "run", objective="Investigate a future capability",
        baseline={"scientific_identity": "sha256:fixed-scientific-environment"},
    )
    monkeypatch.setattr(
        "chem_coworker.scientific_workspace.workspace.verify_baseline", lambda *_, **__: None,
    )
    return ScientificWorkspace(store.root, operations=TrialOperations(store))


def test_new_capability_records_contract_and_verified_links(workspace: ScientificWorkspace) -> None:
    source = workspace.store.append("call", {"result": "observed data"})
    entry = next(item for item in workspace.operations.catalog() if item["name"] == "inspect_trial")
    assert entry["contract_version"] == "2"
    assert entry["evidence_arguments"] == ["observations"]
    event = workspace.run("inspect_trial", {"observations": [source.artifact_ref]})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["operation_contract_version"] == "2"
    assert recorded["scientific_identity"] == "sha256:fixed-scientific-environment"
    assert event.evidence_refs == (source.artifact_ref,)
    assert recorded["execution_status"] == "completed"
    assert recorded["result"]["status"] == "unknown"
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


@pytest.mark.parametrize("status", ["timed_out", "error", "cancelled"])
def test_new_capability_retains_partial_execution_outcome(
    workspace: ScientificWorkspace, status: str,
) -> None:
    event = workspace.run("inspect_trial", {"observations": [], "execution_status": status})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == status
    assert recorded["result"]["status"] == "unknown"
    assert recorded["error"]["type"] == "OperationIncomplete"
    with pytest.raises(ValueError, match="Only completed"):
        workspace.replay(event.artifact_ref)


@pytest.mark.parametrize("result", [
    {"status": "unknown"}, {"execution_status": "running"},
    {"execution_status": None}, {"execution_status": []}, "unfinished response",
])
def test_malformed_declared_status_cannot_become_completed_evidence(
    workspace: ScientificWorkspace, monkeypatch: pytest.MonkeyPatch, result: Any,
) -> None:
    monkeypatch.setattr(workspace.operations, "invoke", lambda *_, **__: result)
    event = workspace.run("inspect_trial", {"observations": []})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == "error"
    assert recorded["result"] == result
    assert "invalid execution_status" in recorded["error"]["message"]
    with pytest.raises(ValueError, match="Only completed"):
        workspace.replay(event.artifact_ref)


def test_replay_rejects_changed_api_before_invocation(workspace: ScientificWorkspace) -> None:
    source = workspace.run("inspect_trial", {"observations": []})

    class UpdatedOperations(TrialOperations):
        DEFINITIONS = tuple(
            replace(item, contract_version="3") if item.name == "inspect_trial" else item
            for item in TrialOperations.DEFINITIONS
        )

        def inspect_trial(self, **arguments: Any) -> Any:
            raise AssertionError("Changed contracts must never replay the old call")

    updated = ScientificWorkspace(workspace.store.root, operations=UpdatedOperations(workspace.store))
    with pytest.raises(ValueError, match="Operation contract changed"):
        updated.replay(source.artifact_ref)
    assert updated.store.read_artifact(source.artifact_ref)["operation_contract_version"] == "2"


def test_historical_calls_without_contract_version_remain_replayable(workspace: ScientificWorkspace) -> None:
    source = workspace.store.append("call", {
        "operation": "inspect_trial", "arguments": {"observations": []},
        "execution_status": "completed",
        "result": {"observations": [], "status": "unknown", "execution_status": "completed"},
    })
    replay = workspace.replay(source.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_unverified_argument_links_are_retained_in_request_only(workspace: ScientificWorkspace) -> None:
    event = workspace.run("inspect_trial", {"observations": ["sha256:missing"]})
    assert event.evidence_refs == ()
    assert workspace.store.read_artifact(event.artifact_ref)["arguments"]["observations"] == ["sha256:missing"]


@pytest.mark.parametrize("name", ["__import__", "catalog", "_path", "not_registered"])
def test_core_rejects_unregistered_or_private_invocation(workspace: ScientificWorkspace, name: str) -> None:
    event = workspace.run(name, {})
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["execution_status"] == "error"
    assert "Unknown scientific operation" in result["error"]["message"]


def test_registry_rejects_duplicate_declarations(workspace: ScientificWorkspace) -> None:
    class DuplicateOperations(TrialOperations):
        DEFINITIONS = TrialOperations.DEFINITIONS + (OperationDefinition("inspect_trial"),)

    with pytest.raises(ValueError, match="unique"):
        DuplicateOperations(workspace.store)


def test_provider_mutation_does_not_change_saved_inputs(workspace: ScientificWorkspace) -> None:
    class MutatingOperations(TrialOperations):
        def invoke(self, operation: str, arguments: Mapping[str, Any]) -> Any:
            arguments["observations"].append("provider-generated")
            return {"observations": arguments["observations"]}

    recorder = ScientificWorkspace(workspace.store.root, operations=MutatingOperations(workspace.store))
    arguments: dict[str, Any] = {"observations": []}
    event = recorder.run("inspect_trial", arguments)
    recorded = recorder.store.read_artifact(event.artifact_ref)
    assert arguments == recorded["arguments"] == {"observations": []}
    assert recorded["result"]["observations"] == ["provider-generated"]
