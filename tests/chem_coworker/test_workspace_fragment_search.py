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


def test_query_preview_and_saved_source_investigation_are_recorded_and_replayable(workspace):
    preview = workspace.run("propose_fragment_queries", {"target_smiles": "CCOC", "query": "COC"})
    assert workspace.store.read_artifact(preview.artifact_ref)["execution_status"] == "completed"
    search = workspace.run("search_fragment_precedents", {"query": "COC", "target_smiles": "CCOC"})
    inspection = workspace.run("investigate_fragment_precedent", {
        "source_ref": search.artifact_ref, "observation_id": "test-observation",
    })
    payload = workspace.store.read_artifact(inspection.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    assert payload["result"]["source_ref"] == search.artifact_ref
    assert payload["result"]["source"]["procedures"]
    assert payload["result"]["transfer"]["status"] == "source_compilation_rejected"
    summary = workspace.call_summary(inspection)["result_summary"]
    assert summary["source"]["observation_id"] == "test-observation"
    for event in (preview, inspection):
        replay = workspace.replay(event.artifact_ref)
        assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_investigation_rejects_forged_source_or_non_search_call(workspace):
    preview = workspace.run("propose_fragment_queries", {"target_smiles": "CCOC", "query": "COC"})
    invalid = workspace.run("investigate_fragment_precedent", {"source_ref": preview.artifact_ref, "observation_id": "test-observation"})
    assert workspace.store.read_artifact(invalid.artifact_ref)["execution_status"] == "error"


def test_automatic_discovery_can_feed_the_recorded_investigation_tool(workspace):
    search = workspace.run("find_synthesis_precedents", {"target_smiles": "COC"})
    payload = workspace.store.read_artifact(search.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    assert payload["result"]["hits"][0]["observation_id"] == "test-observation"
    inspection = workspace.run("investigate_fragment_precedent", {
        "source_ref": search.artifact_ref, "observation_id": "test-observation",
    })
    result = workspace.store.read_artifact(inspection.artifact_ref)
    assert result["execution_status"] == "completed", result
    assert result["result"]["source"]["procedures"]
    assert search.artifact_ref in inspection.evidence_refs
    replay = workspace.replay(inspection.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_search_validates_and_records_selected_preview(workspace):
    preview = workspace.run("propose_fragment_queries", {"target_smiles": "CCOC", "query": "COC"})
    variant = workspace.store.read_artifact(preview.artifact_ref)["result"]["variants"][0]
    request = {"target_smiles": "CCOC", **{key: variant[key] for key in ("query", "query_format", "topology")},
               "query_variant_ref": preview.artifact_ref, "query_variant_id": variant["variant_id"]}
    search = workspace.run("search_fragment_precedents", request)
    assert workspace.store.read_artifact(search.artifact_ref)["execution_status"] == "completed"
    assert preview.artifact_ref in search.evidence_refs
    bad = workspace.run("search_fragment_precedents", {**request, "query": "CO"})
    payload = workspace.store.read_artifact(bad.artifact_ref)
    assert payload["execution_status"] == "error"
    assert "differs from" in payload["error"]["message"]
    search = workspace.run("search_fragment_precedents", {"query": "COC", "target_smiles": "CCOC"})
    invalid = workspace.run("investigate_fragment_precedent", {"source_ref": search.artifact_ref, "observation_id": "forged"})
    assert workspace.store.read_artifact(invalid.artifact_ref)["execution_status"] == "error"


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


def test_console_reuses_library_preserves_evidence_and_records_inspection(workspace):
    from chem_coworker.scientific_workspace.adapters.fragment_search import FragmentWorker
    from chem_coworker.scientific_workspace.adapters.operations import _fragment_replay_result
    from chem_coworker.scientific_workspace.console import ScientificConsole

    arguments = {"query": "COC", "target_smiles": "CCOC"}
    cold = workspace.run("search_fragment_precedents", arguments)
    cold_result = workspace.store.read_artifact(cold.artifact_ref)["result"]
    with FragmentWorker() as worker:
        workspace.operations.fragment_worker = worker
        console = ScientificConsole(workspace)
        first = console.dispatch({"operation": "search_fragment_precedents", "arguments": arguments,
                                  "reason": "Find ether construction"})
        first_ref = first["event"]["artifact_ref"]
        preview = console.dispatch({"operation": "propose_fragment_queries", "arguments": arguments,
                                    "parent_ref": first_ref, "reason": "Check the saved query semantics"})
        second = console.dispatch({"action": "search_prepared", "query_ref": preview["event"]["artifact_ref"],
                                   "reason": "Replay the exact saved parent query"})
        for response in (first, second):
            result = workspace.store.read_artifact(response["event"]["artifact_ref"])["result"]
            assert _fragment_replay_result(result) == _fragment_replay_result(cold_result)
            assert result["execution"]["candidate_set_reused"] is False
        assert first["search_execution"]["library_reused"] is False
        assert second["search_execution"]["library_reused"] is True
        assert second["search_execution"]["session_library_loads"] == 1
        assert first["search_execution"]["session_id"] == second["search_execution"]["session_id"]
        inspected = console.dispatch({"action": "inspect", "source_ref": first_ref,
                                      "path": ["result", "hits", 0, "matches"]})
        assert workspace.store.read_artifact(inspected["inspection_ref"])["review_status"] == "opened_not_adjudicated"
        process = worker.process
    assert process.poll() is not None


def test_warm_worker_deadline_kills_and_next_call_restarts(workspace, monkeypatch):
    import chem_coworker.scientific_workspace.adapters.fragment_search as adapter

    original = adapter._session_worker_command
    with adapter.FragmentWorker() as worker:
        workspace.operations.fragment_worker = worker
        monkeypatch.setattr(adapter, "_session_worker_command", lambda path: [sys.executable, "-c", "import time; time.sleep(60)"])
        started = monotonic()
        event = workspace.run("search_fragment_precedents", {"query": "CO", "timeout_seconds": 1})
        assert monotonic() - started < 8
        payload = workspace.store.read_artifact(event.artifact_ref)
        assert payload["execution_status"] == "timed_out"
        assert payload["result"]["counts"]["products"]["precision"] == "unknown"
        assert worker.process is None
        monkeypatch.setattr(adapter, "_session_worker_command", original)
        event = workspace.run("search_fragment_precedents", {"query": "COC"})
        result = workspace.store.read_artifact(event.artifact_ref)
        assert result["execution_status"] == "completed"
        assert result["result"]["execution"]["library_reused"] is False


def test_console_rejects_invalid_saved_query_and_requires_reason(workspace):
    from chem_coworker.scientific_workspace.console import ScientificConsole

    console = ScientificConsole(workspace)
    with pytest.raises(ValueError, match="chemical search question"):
        console.dispatch({"operation": "search_fragment_precedents", "arguments": {"query": "CO"}})
    bad = workspace.run("propose_fragment_queries", {"query": "c1ccccc1", "target_smiles": "COC"})
    with pytest.raises(ValueError, match="completed query preview"):
        console.dispatch({"action": "search_prepared", "query_ref": bad.artifact_ref, "reason": "Must reject"})


def test_lazy_session_validates_target_before_opening_index(tmp_path):
    from condition_recommender.fragment_search import FragmentSearchSession, search_fragment_precedents

    missing = tmp_path / "missing.sqlite"
    session = FragmentSearchSession(missing)
    with pytest.raises(ValueError, match="does not match target_smiles"):
        search_fragment_precedents(missing, "c1ccccc1", target_smiles="COC", _session=session)
    assert session.connection is None
    session.__exit__(None, None, None)


def test_session_cannot_search_another_index(workspace, tmp_path):
    from condition_recommender.fragment_search import FragmentSearchSession, search_fragment_precedents

    path = workspace.store.manifest["baseline"]["artifacts"]["fragment_index"]["path"]
    with FragmentSearchSession(path) as session:
        with pytest.raises(ValueError, match="different index"):
            search_fragment_precedents(tmp_path / "another.sqlite", "CO", _session=session)
    assert session.connection is None
    assert session.library is None
    assert session.candidates == {}
