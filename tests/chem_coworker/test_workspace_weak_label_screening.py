"""Recorded weak-label screening, dataset identity and agent conversation integration."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from fastapi.testclient import TestClient
import pytest

from app.web_api.main import create_app
from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.runtime.agent_runtime import AgentResult
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity,
    code_manifest,
    environment_versions,
)
from chem_coworker.scientific_workspace.views.call_summaries import summarize_call
from chem_coworker.scientific_workspace.views.brief_summaries import summarize_call_brief
from chem_coworker.scientific_workspace.runtime.conversation import ConversationService
from condition_recommender import (
    generate_weak_label_screening_array,
    weak_label_recipe_catalog_path,
)
from tests.chem_coworker.test_scientific_conversation import finish
from tests.condition_recommender.test_weak_label_recommendation import QUERY, _dataset


ROOT = Path(__file__).resolve().parents[2]
OPERATION = "generate_weak_label_screening_array"


@pytest.fixture
def workspace(tmp_path: Path) -> ScientificWorkspace:
    records = _dataset(tmp_path)
    baseline = {
        "repository": str(ROOT), "code_files": code_manifest(ROOT),
        "environment": environment_versions(),
        "artifacts": {
            "weak_label_records": artifact_identity(records),
            "weak_label_recipe_catalog": artifact_identity(weak_label_recipe_catalog_path(records)),
        },
    }
    InvestigationStore.create(tmp_path / "investigation", objective="Screen conditions", baseline=baseline)
    return ScientificWorkspace(tmp_path / "investigation")


def test_screening_preserves_domain_output_without_structural_indexes_and_replays(
    workspace: ScientificWorkspace,
) -> None:
    assert OPERATION in {item["name"] for item in workspace.operations.catalog()}
    event = workspace.run(OPERATION, {"reaction_smiles": QUERY, "array_size": 96})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == "completed", recorded
    records = workspace.store.manifest["baseline"]["artifacts"]["weak_label_records"]["path"]
    expected = generate_weak_label_screening_array(QUERY, records_path=records, array_size=96)
    assert recorded["result"] == json.loads(json.dumps(expected.to_dict()))
    assert recorded["result"]["valid"]
    assert len(recorded["result"]["recommendations"]) == 2  # Never pad a short pool.
    assert "WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED" in recorded["result"]["warnings"]
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"]
    capability = workspace.capabilities()["weak_label_screening"]
    assert capability["requires_condition_index"] is False
    assert capability["source_structures_verified"] is False


@pytest.mark.parametrize("size", [0, -1, 251, True, 24.0, "96"])
def test_invalid_sizes_are_recorded_errors(workspace: ScientificWorkspace, size: Any) -> None:
    event = workspace.run(OPERATION, {"reaction_smiles": QUERY, "array_size": size})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == "error"
    assert "array_size must be an integer" in recorded["error"]["message"]


def test_maximum_size_and_default_are_supported(workspace: ScientificWorkspace) -> None:
    for arguments in ({"reaction_smiles": QUERY}, {"reaction_smiles": QUERY, "array_size": 250}):
        event = workspace.run(OPERATION, arguments)
        assert workspace.store.read_artifact(event.artifact_ref)["result"]["valid"]


@pytest.mark.parametrize("reaction,hint,valid", [
    (QUERY, "SNAr", True),
    (QUERY, "Suzuki-Miyaura", False),
    ("CCO>>CCO", None, False),
    ("[CH3:1][Br:2].[NH3:1]>>[CH3:1][NH2:1]", None, False),
])
def test_query_evidence_and_conflicting_labels_remain_domain_results(
    workspace: ScientificWorkspace, reaction: str, hint: str | None, valid: bool,
) -> None:
    event = workspace.run(OPERATION, {"reaction_smiles": reaction, "source_reaction_type_hint": hint})
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["execution_status"] == "completed"
    assert result["result"]["valid"] is valid
    if valid:
        assert result["result"]["source_reaction_type_candidates"] == ["SNAr"]
    else:
        assert result["result"]["error"]
        assert result["result"]["recommendations"] == []


@pytest.mark.parametrize("missing", ["weak_label_records", "weak_label_recipe_catalog"])
def test_missing_inputs_are_explicit_and_never_use_an_unpinned_default(
    workspace: ScientificWorkspace, missing: str,
) -> None:
    del workspace.store.manifest["baseline"]["artifacts"][missing]
    event = workspace.run(OPERATION, {"reaction_smiles": QUERY})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == "error"
    assert missing in recorded["error"]["message"]


@pytest.mark.parametrize("changed", ["weak_label_records", "weak_label_recipe_catalog"])
def test_changed_inputs_block_new_calls(workspace: ScientificWorkspace, changed: str) -> None:
    path = Path(workspace.store.manifest["baseline"]["artifacts"][changed]["path"])
    with path.open("ab") as handle:
        handle.write(b"changed")
    event = workspace.run(OPERATION, {"reaction_smiles": QUERY})
    recorded = workspace.store.read_artifact(event.artifact_ref)
    assert recorded["execution_status"] == "error"
    assert "Baseline artifact changed" in recorded["error"]["message"]


def test_creation_pins_sibling_catalog_and_reports_missing_catalog(tmp_path: Path) -> None:
    records = _dataset(tmp_path)
    catalog = weak_label_recipe_catalog_path(records)
    workspace = ScientificWorkspace.create(
        tmp_path / "created", objective="Pin both screening inputs", repository=ROOT,
        artifacts={"weak_label_records": records},
    )
    inputs = workspace.store.manifest["baseline"]["artifacts"]
    assert inputs["weak_label_recipe_catalog"] == artifact_identity(catalog)
    catalog.unlink()
    missing = ScientificWorkspace.create(
        tmp_path / "missing", objective="Report missing screening input", repository=ROOT,
        artifacts={"weak_label_records": records},
    )
    assert missing.store.manifest["baseline"]["artifacts"]["weak_label_recipe_catalog"]["status"] == "missing"
    event = missing.run(OPERATION, {"reaction_smiles": QUERY})
    assert "weak_label_recipe_catalog" in missing.store.read_artifact(event.artifact_ref)["error"]["message"]


def test_creation_rejects_an_unrelated_catalog(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="catalog beside weak_label_records"):
        ScientificWorkspace.create(
            tmp_path / "bad", objective="Prevent unpinned catalog use", repository=ROOT,
            artifacts={"weak_label_records": _dataset(tmp_path),
                       "weak_label_recipe_catalog": tmp_path / "other.gz"},
        )


def test_summaries_count_full_array_and_retain_weak_label_provenance(workspace: ScientificWorkspace) -> None:
    event = workspace.run(OPERATION, {"reaction_smiles": QUERY})
    payload = workspace.store.read_artifact(event.artifact_ref)
    payload["result"]["recommendations"] *= 48  # Exercise a large saved result projection.
    for project in (summarize_call, summarize_call_brief):
        summary = project(payload)["result_summary"]
        assert summary["recommendation_mode"] == "weak_label_screening"
        assert summary["recommendation_count"] == 96
        assert len(summary["recommendations"]) == 3
        assert "WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED" in summary["warnings"]
        assert summary["recommendations"][0]["source_row_numbers"]
        assert summary["recommendations"][0]["cautions"]


def test_agent_web_conversation_records_screening_and_boots(tmp_path: Path) -> None:
    class ScreeningRuntime:
        """Substitute the model while exercising real conversation and chemistry code."""

        def describe(self) -> dict[str, Any]:
            return {"runtime": "test_double"}

        def run(self, **kwargs: Any) -> AgentResult:
            assert OPERATION in kwargs["prompt"]
            scientific = ScientificWorkspace(kwargs["workspace"])
            assert "w.task_guide('conditions')" in kwargs["prompt"]
            assert "conditions_screening" in scientific.task_guide("conditions")["text"]
            guide = scientific.task_guide("conditions_screening")
            assert "source row" in guide["text"]
            frozen = scientific.store.manifest["baseline"]["learning_context"]
            assert guide == frozen["guides"]["conditions_screening"]
            event = scientific.run(OPERATION, {"reaction_smiles": QUERY, "array_size": 96})
            result = scientific.store.read_artifact(event.artifact_ref)
            assert result["result"]["valid"]
            return AgentResult({
                "schema_version": "scientific_answer.v2",
                "sources": [], "molecules": [], "target_molecule_ids": [],
                "steps": [], "routes": [], "claims": [],
                "answer_markdown": "Two weak-label screening recipes; source structures are unverified.",
                "evidence_refs": [event.artifact_ref], "uncertainties": ["Not experimentally validated"],
                "needs_user_input": False,
            }, "screening-thread", {})

    records = _dataset(tmp_path)
    service = ConversationService(tmp_path / "conversations", runtime=ScreeningRuntime(),
                                  artifacts={"weak_label_records": records})
    try:
        client = TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False),
                            base_url="http://127.0.0.1")
        assert client.get("/scientific").status_code == 200
        identity = service.submit("Generate 96 diverse screening recipes")["conversation_id"]
        turn = finish(service, identity)
        assert turn["status"] == "completed", turn
        evidence = service.artifact(identity, turn["answer"]["evidence_refs"][0])
        assert evidence["operation"] == OPERATION
        assert evidence["result"]["recommendation_mode"] == "weak_label_screening"
    finally:
        service.close()


def test_default_agent_configuration_includes_weak_label_records() -> None:
    config = json.loads((ROOT / "examples/ai_native/artifacts.local.example.json").read_text("utf-8"))
    assert config["weak_label_records"] == "datasets/weak_label/v2.1_cleaned.csv"
