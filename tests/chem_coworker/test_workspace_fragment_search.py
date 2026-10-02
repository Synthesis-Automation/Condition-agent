"""Recorded fragment search, bounded inspection, hard deadlines, and replay."""

from dataclasses import asdict
import json
from pathlib import Path
import sys
from time import monotonic

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity,
    code_manifest,
    environment_versions,
)
from condition_recommender.fragment_index import build_fragment_index
from reactive_taxonomy import featurize_reaction

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def workspace(tmp_path):
    reaction = "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]"
    source = tmp_path / "records.jsonl"
    source.write_text(json.dumps({"observation_id": "test-observation", "reaction_id": "reaction-1",
                                  "reference_id": "patent-1", "reaction_smiles": reaction,
                                  "admission_tier": "review", "admission_reasons": ["missing_conditions"],
                                  "reaction_observation": asdict(featurize_reaction(reaction).observation)}), "utf-8")
    catalog = tmp_path / "procedures.jsonl"
    catalog.write_text(json.dumps({"observation_id": "test-observation", "reaction_id": "reaction-1",
                                   "text": "Reported source text. " * 100}), "utf-8")
    index = tmp_path / "fragment.sqlite"
    build_fragment_index(source, index, procedure_catalog=catalog)
    baseline = {"repository": str(ROOT), "code_files": code_manifest(ROOT),
                "environment": environment_versions(), "artifacts": {"fragment_index": artifact_identity(index)}}
    InvestigationStore.create(tmp_path / "workspace", objective="Inspect an unfamiliar core", baseline=baseline)
    return ScientificWorkspace(tmp_path / "workspace")


def test_real_worker_compact_summary_evidence_inspection_and_replay(workspace):
    event = workspace.run("search_fragment_precedents", {"query": "COC", "target_smiles": "CCOC"})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    result = payload["result"]
    assert result["search_status"] == "complete"
    assert result["hits"][0]["admission_tier"] == "review"
    summary = workspace.call_summary(event)["result_summary"]
    assert summary["search_status"] == "complete"
    assert summary["hits"][0]["product_smiles"] == "COC"
    assert summary["hits"][0]["procedure_availability"] == "linked"
    assert summary["target_validation"]["matches_target"] is True
    assert summary["target_validation"]["target_smiles"] == "CCOC"
    path = ("result", "hits", 0, "procedures", 0, "record", "text", "chunks")
    page = workspace.inspect_artifact(event.artifact_ref, path, offset=1, limit=2)
    assert page["preview"][0]["start"] == 300
    assert page["page"]["next_offset"] == 3
    stages = Path(result["execution"]["diagnostics"]["stages.jsonl"]).read_text("utf-8")
    assert "matching_and_evidence_finished" in stages
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_invalid_query_is_error_not_zero_hits(workspace):
    event = workspace.run("search_fragment_precedents", {"query": "C.O"})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "error"
    assert payload["error"]["type"] == "ValueError"
    assert "connected" in payload["error"]["message"]


def test_target_mismatch_is_recorded_as_error_not_absence(workspace):
    event = workspace.run("search_fragment_precedents", {
        "query": "O=C1CCCC2OC3CCC(C3)N12",
        "target_smiles": "O=C1c2ccccc2C[C@H]3O[C@@H](C4)CC[C@@H]4N13",
    })
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "error"
    assert "does not match target_smiles" in payload["error"]["message"]
    assert "counts" not in payload["result"]


def test_worker_deadline_stops_child_and_saves_diagnostics(workspace, monkeypatch):
    import chem_coworker.scientific_workspace.adapters.fragment_search as adapter
    monkeypatch.setattr(adapter, "_worker_command", lambda path: [sys.executable, "-c", "import time; time.sleep(60)"])
    started = monotonic()
    event = workspace.run("search_fragment_precedents", {"query": "CO", "timeout_seconds": 1})
    assert monotonic() - started < 8
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "timed_out"
    assert payload["result"]["search_status"] == "partial"
    assert payload["result"]["counts"]["products"]["precision"] == "unknown"
    assert "stderr.log" in payload["result"]["execution"]["diagnostics"]
    with pytest.raises(ValueError, match="completed"):
        workspace.replay(event.artifact_ref)


def test_changed_index_is_rejected_before_search(workspace):
    index = Path(workspace.store.manifest["baseline"]["artifacts"]["fragment_index"]["path"])
    with index.open("ab") as handle:
        handle.write(b"modified")
    event = workspace.run("search_fragment_precedents", {"query": "COC"})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "error"
    assert "artifact changed" in payload["error"]["message"]


def test_suggestion_is_recorded_replayable_and_requires_no_index(tmp_path, monkeypatch):
    import condition_recommender.fragment_search as search

    monkeypatch.setattr(search, "search_fragment_precedents", lambda *a, **kw: pytest.fail("Unexpected search"))
    baseline = {"repository": str(ROOT), "code_files": code_manifest(ROOT),
                "environment": environment_versions(), "artifacts": {}}
    InvestigationStore.create(tmp_path / "workspace", objective="Choose an unfamiliar core", baseline=baseline)
    workspace = ScientificWorkspace(tmp_path / "workspace")
    event = workspace.run("suggest_search_fragments", {"target_smiles": "CC(=O)c1ccc2c(c1)COc1ccccc1-2"})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    result = payload["result"]
    assert result["candidates"][0]["query"] == "c1ccc2c(c1)COc1ccccc1-2"
    assert result["candidates"][0]["boundaries"]
    summary = workspace.call_summary(event)["result_summary"]
    assert summary["candidates"][0]["query"] == result["candidates"][0]["query"]
    assert summary["definition_version"] == "search_fragments.v1@1.0"
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True
    invalid = workspace.run("suggest_search_fragments", {"target_smiles": "C.O"})
    assert workspace.store.read_artifact(invalid.artifact_ref)["execution_status"] == "error"
