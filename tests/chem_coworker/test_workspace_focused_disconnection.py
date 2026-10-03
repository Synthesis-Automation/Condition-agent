"""Focused single-step tools retain target binding, evidence and replay."""

from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity, code_manifest, environment_versions,
)
from core_retrosynthesis.generic_library import save_generic_library
from tests.core_retrosynthesis_tests.test_focused_disconnection import library


@pytest.fixture
def workspace(tmp_path, library):
    path = tmp_path / "library.json.gz"
    save_generic_library(library, path)
    root = Path(__file__).resolve().parents[2]
    baseline = {"repository": str(root), "code_files": code_manifest(root),
                "environment": environment_versions(),
                "artifacts": {"retro_library": artifact_identity(path)}}
    InvestigationStore.create(tmp_path / "workspace", objective="Focus one bond", baseline=baseline)
    return ScientificWorkspace(tmp_path / "workspace")


def test_focused_tool_record_summary_and_replay(workspace):
    event = workspace.run("disconnect_target", {
        "target_smiles": "CNCC", "required_disconnection_bond": [2, 3],
        "focus_target_smiles": "CCNC", "top_k": 1, "max_candidates_to_validate": 1,
    })
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    result = payload["result"]
    assert result["valid"] and result["bond_focus"]["atom_ids"] == [2, 3]
    candidate = result["strategies"][0]["representative"]
    assert candidate["precursor_smiles"] == "CBr.CCN"
    assert candidate["bond_focus_check"]["status"] == "verified"
    assert result["search_diagnostics"]["validation_attempt_count"] == 1
    assert payload["operation_contract_version"] == "2"
    for detailed in (True, False):
        summary = workspace.call_summary(event, detailed=detailed)["result_summary"]
        assert summary["bond_focus"]["atom_ids"] == [2, 3]
        assert summary["search_diagnostics"]
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"]


def test_stale_focus_and_no_match_have_distinct_outcomes(workspace):
    invalid = workspace.run("disconnect_target", {
        "target_smiles": "CCNC", "required_disconnection_bond": [1, 2],
        "focus_target_smiles": "CCO",
    })
    invalid_result = workspace.store.read_artifact(invalid.artifact_ref)["result"]
    assert not invalid_result["valid"] and "identify" in invalid_result["error"]
    absent = workspace.run("disconnect_target", {
        "target_smiles": "CCNC", "required_disconnection_bond": [0, 1],
        "focus_target_smiles": "CCNC",
    })
    result = workspace.store.read_artifact(absent.artifact_ref)["result"]
    assert result["valid"] and not result["strategies"]
    assert result["bond_focus"]["atom_ids"] == [0, 1]
    assert result["search_diagnostics"]["validation_attempt_count"] == 0
    assert "bounded" in result["warnings"][0]
