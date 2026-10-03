"""Starting-material discovery, summaries, evidence recording and replay."""

import json
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity, code_manifest, environment_versions,
)
from chem_coworker.scientific_workspace.views.brief_summaries import summarize_call_brief
from condition_recommender.fragment_index import build_fragment_index


@pytest.fixture
def workspace(tmp_path):
    source = tmp_path / "source.jsonl"
    source.write_text(json.dumps({"observation_id": "obs-1", "reaction_id": "rxn-1",
                                  "reference_id": "ref-1", "reaction_smiles": "CC=O>>CCO"}), "utf-8")
    index = tmp_path / "index.sqlite"
    build_fragment_index(source, index)
    root = Path(__file__).resolve().parents[2]
    baseline = {"repository": str(root), "code_files": code_manifest(root),
                "environment": environment_versions(),
                "artifacts": {"fragment_index": artifact_identity(index)}}
    InvestigationStore.create(tmp_path / "workspace", objective="Assess route leaves", baseline=baseline)
    return ScientificWorkspace(tmp_path / "workspace")


@pytest.mark.parametrize("options,reason", [
    ({}, "registry_match"),
    ({"allow_registry_stop": False}, "exact_literature_match"),
    ({"allow_registry_stop": False, "allow_literature_stop": False}, "low_molecular_weight"),
    ({"unavailable_starting_materials": ["OCC"]}, "explicitly_unavailable"),
])
def test_recorded_decisions_summary_and_replay(workspace, options, reason):
    catalog = {item["name"]: item for item in workspace.operations.catalog()}
    assert catalog["assess_starting_material"]["required_artifacts"] == []
    assert workspace.capabilities()["starting_material_assessment"]["availability_verified"] is False
    event = workspace.run("assess_starting_material", {"smiles": "CCO", **options})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    assert payload["result"]["stop_reason"] == reason
    detailed = workspace.call_summary(event)["result_summary"]
    brief = summarize_call_brief(payload)["result_summary"]
    for summary in (detailed, brief):
        assert summary["stop_reason"] == reason
        assert summary["availability"] == "unknown"
        assert summary["warnings"]
        assert "stop_expansion" in summary
        if reason == "exact_literature_match":
            assert summary["literature"]["precedents"][0]["observation_id"] == "obs-1"
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_invalid_material_remains_recorded_error(workspace):
    event = workspace.run("assess_starting_material", {"smiles": "not_smiles"})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "error"
