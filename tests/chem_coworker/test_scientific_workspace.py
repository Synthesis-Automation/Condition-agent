"""Scientific workspace provenance, recovery, and canonical-operation regressions."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.__main__ import main
from chem_coworker.scientific_workspace.baseline import (
    artifact_identity, code_manifest, environment_versions, verify_baseline,
)
from chem_coworker.scientific_workspace.store import canonical_bytes
from reactive_taxonomy import featurize_reaction


ROOT = Path(__file__).resolve().parents[2]


def test_create_captures_real_code_versions_registry_and_missing_inputs(tmp_path: Path) -> None:
    workspace = ScientificWorkspace.create(
        tmp_path / "created", objective="Check preparation", repository=ROOT,
        artifacts={"condition_index": tmp_path / "missing.sqlite"},
    )
    baseline = workspace.store.manifest["baseline"]
    assert baseline["git_revision"]
    assert baseline["artifacts"]["condition_index"]["status"] == "missing"
    assert "registry_validation" in baseline
    assert baseline["taxonomy_versions"]
    assert baseline["registry_versions"]
    verify_baseline(baseline)


@pytest.fixture
def workspace(tmp_path: Path) -> ScientificWorkspace:
    baseline = {
        "repository": str(ROOT), "code_files": code_manifest(ROOT),
        "environment": environment_versions(), "artifacts": {},
        "validation_status": "development_snapshot_not_release_validated",
    }
    InvestigationStore.create(tmp_path / "investigation", objective="Investigate chemistry",
                              baseline=baseline, constraints=("Preserve uncertainty",))
    return ScientificWorkspace(tmp_path / "investigation")


def test_reopen_preserves_evidence_decisions_and_lifecycle(workspace: ScientificWorkspace) -> None:
    event = workspace.run("analyze_reaction", {"reaction_smiles": "CCBr.N>>CCN"})
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["execution_status"] == "completed"
    workspace.store.note("hypothesis", "Transfer remains uncertain", evidence_refs=(event.artifact_ref,))
    workspace.store.set_status("insufficient_evidence", "Need a procedure")
    reopened = ScientificWorkspace(workspace.store.root)
    summary = reopened.store.summary()
    assert summary["status"] == "insufficient_evidence"
    assert summary["constraints"] == ["Preserve uncertainty"]
    assert summary["notes"][0]["origin"] == "agent_authored"
    assert summary["notes"][0]["review_status"] == "unreviewed"
    assert len(summary["events"]) == 3
    stopped = reopened.run("analyze_molecule", {"smiles": "CCO"})
    assert reopened.store.read_artifact(stopped.artifact_ref)["execution_status"] == "error"
    reopened.store.set_status("active", "Procedure now available")
    replay = reopened.replay(event.artifact_ref)
    assert reopened.store.read_artifact(replay.artifact_ref)["matches"] is True


@pytest.mark.parametrize("reaction", [
    "CCBr.N>>CCN", "N.CCBr>>CCN", "CC>>CC", "not_a_reaction",
    "[CH3:1][Br:2].[NH3:1]>>[CH3:1][NH2:1]",
])
def test_reaction_adapter_preserves_complete_domain_result(
    workspace: ScientificWorkspace, reaction: str,
) -> None:
    event = workspace.run("analyze_reaction", {"reaction_smiles": reaction})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == "completed"
    assert canonical_bytes(recorded["result"]) == canonical_bytes(featurize_reaction(reaction))


def test_failed_calls_are_recorded_and_do_not_destroy_previous_results(workspace: ScientificWorkspace) -> None:
    good = workspace.run("analyze_molecule", {"smiles": "CCO"})
    missing = workspace.run("recommend_conditions", {"reaction_smiles": "CCBr.N>>CCN"})
    unsupported = workspace.run("__import__", {"name": "os"})
    invalid_ref = workspace.run("revise_routes", {"source_ref": "../bad", "intent": {}})
    assert workspace.store.read_artifact(good.artifact_ref)["execution_status"] == "completed"
    for event in (missing, unsupported, invalid_ref):
        assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "error"
    assert "condition_index" in workspace.store.read_artifact(missing.artifact_ref)["error"]["message"]


def test_custom_files_are_snapshotted_without_execution(workspace: ScientificWorkspace, tmp_path: Path) -> None:
    script = tmp_path / "analysis.py"
    script.write_text("raise RuntimeError('must not execute')\n", "utf-8")
    event = workspace.store.attach_file(script, description="Unreviewed comparison")
    script.write_text("changed", "utf-8")
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["content"].startswith("raise RuntimeError")
    assert result["origin"] == "derived_analysis"


def test_recipe_operations_preserve_registry_and_compatibility_contracts(workspace: ScientificWorkspace) -> None:
    from condition_registry import ConditionComponentInput, build_resolved_recipe_from_inputs
    from condition_recommender import assess_reaction_recipe

    components = [{"raw_identifier": "ethanol", "identifier_type": "name", "source_field": "user"}]
    resolved = workspace.run("resolve_recipe", {"components": components})
    result = workspace.store.read_artifact(resolved.artifact_ref)["result"]
    assert canonical_bytes(result) == canonical_bytes(build_resolved_recipe_from_inputs(
        [ConditionComponentInput(**item) for item in components],
    ))
    for reaction in ("CCBr.N>>CCN", "CC>>CC"):
        event = workspace.run("assess_recipe", {"reaction_smiles": reaction, "recipe": result})
        assert canonical_bytes(workspace.store.read_artifact(event.artifact_ref)["result"]) == canonical_bytes(
            assess_reaction_recipe(reaction, result),
        )


def test_interrupt_records_cancelled_call_and_allows_explicit_resume(workspace: ScientificWorkspace, monkeypatch) -> None:
    def interrupted(operation, arguments):
        raise KeyboardInterrupt()

    monkeypatch.setattr(workspace.operations, "invoke", interrupted)
    event = workspace.run("analyze_molecule", {"smiles": "CCO"})
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "cancelled"
    assert workspace.store.summary()["status"] == "cancelled"
    workspace.store.set_status("active", "Resume after interruption")
    assert workspace.store.summary()["status"] == "active"


def test_artifact_tampering_and_invalid_references_are_rejected(workspace: ScientificWorkspace) -> None:
    event = workspace.store.note("question", "What evidence is missing?")
    with pytest.raises(ValueError, match="Invalid artifact"):
        workspace.store.read_artifact("../../outside")
    with pytest.raises(ValueError, match="Invalid artifact"):
        workspace.store.note("decision", "Unsupported reference", evidence_refs=("bad",))
    path = workspace.store.root / "artifacts" / (event.artifact_ref.split(":")[1] + ".json")
    path.write_text("{}", "utf-8")
    with pytest.raises(ValueError, match="checksum"):
        workspace.store.summary()


def test_existing_investigation_is_never_overwritten(workspace: ScientificWorkspace) -> None:
    with pytest.raises(FileExistsError):
        InvestigationStore.create(workspace.store.root, objective="replacement", baseline={})
    assert workspace.store.summary()["objective"] == "Investigate chemistry"


def test_writer_lock_prevents_conflicting_event_appends(workspace: ScientificWorkspace) -> None:
    lock = workspace.store.root / ".writer.lock"
    lock.write_text("another writer", "utf-8")
    with pytest.raises(RuntimeError, match="writer active"):
        workspace.store.note("question", "Concurrent write")
    assert workspace.store.events() == ()


def test_replay_detects_content_drift_even_when_file_stat_is_unchanged(tmp_path: Path) -> None:
    artifact = tmp_path / "data.json"
    artifact.write_text('{"value":1}', "utf-8")
    identity = artifact_identity(artifact)
    baseline = {"repository": str(tmp_path), "code_files": {},
                "environment": environment_versions(), "artifacts": {"data": identity}}
    artifact.write_text('{"value":2}', "utf-8")
    os.utime(artifact, ns=(identity["mtime_ns"], identity["mtime_ns"]))
    verify_baseline(baseline)
    with pytest.raises(ValueError, match="content changed"):
        verify_baseline(baseline, full_hash=True)


def test_changed_scientific_code_blocks_further_calculations(workspace: ScientificWorkspace) -> None:
    workspace.store.manifest["baseline"]["code_files"] = {}
    event = workspace.run("analyze_molecule", {"smiles": "CCO"})
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["execution_status"] == "error"
    assert "code or definitions changed" in result["error"]["message"]


def test_procedure_access_preserves_multiple_observations_and_missingness(
    workspace: ScientificWorkspace, tmp_path: Path,
) -> None:
    catalog = tmp_path / "procedures.jsonl.gz"
    observations = [
        {"reaction_id": "r1", "observation_id": "o1", "procedure_text": "Reported procedure"},
        {"reaction_id": "r1", "observation_id": "o2", "procedure_text": None},
        {"reaction_id": "r2", "observation_id": "o3", "procedure_text": "Other reaction"},
    ]
    with gzip.open(catalog, "wt", encoding="utf-8") as handle:
        for item in observations:
            handle.write(json.dumps(item) + "\n")
    workspace.store.manifest["baseline"]["artifacts"]["procedure_catalog"] = artifact_identity(catalog)
    result = workspace.operations.get_procedures(["r1", "missing"])
    assert result["records"] == observations[:2]
    assert result["missing_reaction_ids"] == ["missing"]


def test_precedent_pagination_preserves_observation_identity(workspace: ScientificWorkspace) -> None:
    @dataclass(frozen=True)
    class Row:
        reaction_id: str
        observation_id: str
        condition_status: str

    rows = (Row("r1", "o1", "resolved"), Row("r1", "o2", "uncertain"))
    index = SimpleNamespace(reaction_ids={"r1": (0, 1)},
                            select=lambda positions: [rows[index] for index in positions],
                            precedent_scope=SimpleNamespace(value="trusted"))
    workspace.operations._recommender = SimpleNamespace(index=index)
    first = workspace.operations.get_precedents(["r1", "r1", "missing"], limit=1)
    second = workspace.operations.get_precedents(["r1"], offset=first["next_offset"], limit=1)
    assert first["total"] == 2
    assert second["records"][0]["observation_id"] == "o2"
    assert second["records"][0]["condition_status"] == "uncertain"
    assert first["missing_reaction_ids"] == ["missing"]
    assert second["next_offset"] is None


def test_cli_returns_nonzero_for_execution_error_and_can_resume(
    workspace: ScientificWorkspace, tmp_path: Path, capsys: pytest.CaptureFixture[str],
) -> None:
    request = tmp_path / "request.json"
    request.write_text('{"reaction_smiles":"CCBr.N>>CCN"}', "utf-8")
    root = str(workspace.store.root)
    assert main(["run", root, "recommend_conditions", "--input", str(request)]) == 1
    output = json.loads(capsys.readouterr().out)
    assert output["execution_status"] == "error"
    assert main(["summary", root]) == 0
    assert len(json.loads(capsys.readouterr().out)["events"]) == 1
