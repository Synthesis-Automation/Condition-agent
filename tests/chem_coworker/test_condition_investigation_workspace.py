"""Investigation joins, proposals and actual custom execution with immutable evidence."""

from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace
from threading import Event
import subprocess

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.baseline import artifact_identity, code_manifest, environment_versions
from condition_registry import ConditionComponentInput, build_resolved_recipe_from_inputs
from tests.condition_recommender.test_condition_investigation import precedent


@pytest.fixture
def workspace(tmp_path: Path) -> ScientificWorkspace:
    repository = Path(__file__).resolve().parents[2]
    InvestigationStore.create(tmp_path / "investigation", objective="Development condition comparison",
                              baseline={"repository": str(repository), "code_files": code_manifest(repository),
                                        "environment": environment_versions(), "artifacts": {}})
    return ScientificWorkspace(tmp_path / "investigation")


def prepare_inspection(workspace: ScientificWorkspace) -> str:
    components = [ConditionComponentInput("ethanol", "user", identifier_type="name")]
    row = precedent(asdict(build_resolved_recipe_from_inputs(components)))
    index = SimpleNamespace(reaction_ids={"r1": (0,)}, select=lambda positions: [row for _ in positions],
                            precedent_scope=SimpleNamespace(value="trusted"))
    workspace.operations._recommender = SimpleNamespace(index=index)
    return workspace.run("inspect_condition_precedents", {"reaction_smiles": "CCBr.N>>CCN", "reaction_ids": ["r1"]}).artifact_ref


def test_procedure_join_never_borrows_another_observations_conditions(workspace: ScientificWorkspace) -> None:
    catalog = workspace.store.root / "procedures.jsonl"
    records = [{"reaction_id": "r1", "observation_id": identity, "procedure_text": text,
                "source": {"locator": "Example " + str(index)}}
               for index, (identity, text) in enumerate((("o1", None), ("o2", "Other procedure"), (None, "Unassigned procedure")))]
    catalog.write_text("\n".join(json.dumps(row) for row in records), "utf-8")
    workspace.store.manifest["baseline"]["artifacts"]["procedure_catalog"] = artifact_identity(catalog)
    reference = prepare_inspection(workspace)
    result = workspace.store.read_artifact(reference)["result"]
    row = result["precedents"][0]
    assert row["procedure_observations"] == records[:1]
    assert row["reaction_level_procedures"] == records[2:]
    assert len(result["procedure_catalog"]["records"]) == 3


def proposal_arguments(reference: str) -> dict:
    return dict(source_ref=reference, observation_id="o1",
                components=[{"raw_identifier": "ethanol", "source_field": "user", "identifier_type": "name"}],
                operating_conditions={"temperature_c": 30.0}, change_reasons={"temperature_c": "Development hypothesis only"},
                evidence_refs=[reference], assumptions=["Source inspection is not support for the invented temperature"],
                risks=["No experimental transfer has been established"])


def test_adaptation_is_separate_proposal_preserving_original_and_canonical_assessment(workspace: ScientificWorkspace) -> None:
    reference = prepare_inspection(workspace)
    before = workspace.store.read_artifact(reference)
    event = workspace.run("propose_condition_adaptation", proposal_arguments(reference))
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["execution_status"] == "completed", result
    proposal = result["result"]
    assert proposal["origin"] == "agent_proposal"
    assert proposal["transfer_status"] == "not_established"
    assert proposal["changes"][0]["basis"] == "proposed"
    assert proposal["original_recipe"]["temperature_c"] is None
    assert proposal["proposed_recipe"]["temperature_c"] == 30.0
    assert proposal["compatibility"]["schema_version"] == "reaction_recipe_assessment.v2"
    assert workspace.store.read_artifact(reference) == before
    assert reference in event.evidence_refs


@pytest.mark.parametrize("field,value", [("observation_id", "absent"), ("change_reasons", {}), ("evidence_refs", []), ("risks", [])])
def test_adaptation_rejects_unattributed_or_uninspectable_changes(workspace: ScientificWorkspace, field, value) -> None:
    reference = prepare_inspection(workspace)
    arguments = proposal_arguments(reference)
    arguments[field] = value
    event = workspace.run("propose_condition_adaptation", arguments)
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "error"


def test_actual_custom_execution_records_inputs_and_is_usable_as_computed_evidence(workspace: ScientificWorkspace) -> None:
    script = workspace.store.root / "sum.py"
    script.write_text("import json,sys\nfrom pathlib import Path\nx=json.loads(Path(sys.argv[1]).read_text())\nPath(sys.argv[2]).write_text(json.dumps({'sum':sum(x['parameters']['numbers'])}))\n", "utf-8")
    event = workspace.run_python("sum.py", {"numbers": [2, 3]})
    result = workspace.store.read_artifact(event.artifact_ref)
    assert result["execution_status"] == "completed", result
    assert result["result"] == {"sum": 5}
    assert result["script_sha256"] and result["input_sha256"]
    answer = ScientificAnswer.model_validate({
        "schema_version": "scientific_answer.v2", "answer_markdown": "A recorded calculation", "evidence_refs": [],
        "uncertainties": [], "needs_user_input": False, "molecules": [], "target_molecule_ids": [], "steps": [], "routes": [],
        "sources": [{"id": "calc", "kind": "local_artifact", "title": "Custom calculation", "artifact_ref": event.artifact_ref, "url": None, "locator": "result.sum"}],
        "claims": [{"text": "2+3=5", "basis": "computed", "source_ids": ["calc"], "limitations": []}],
    })
    assert validate_answer_evidence(answer, workspace.store) == [event.artifact_ref]
    with pytest.raises(ValueError, match="completed scientific calls"):
        workspace.replay(event.artifact_ref)
    script.write_text("raise RuntimeError('changed original')", "utf-8")
    assert workspace.store.read_artifact(event.artifact_ref)["result"] == {"sum": 5}


@pytest.mark.parametrize("code,status", [
    ("raise RuntimeError('expected failure')", "error"),
    ("import time; time.sleep(10)", "timed_out"),
    ("import sys; from pathlib import Path; Path(sys.argv[2]).write_text('invalid JSON')", "error"),
    ("import sys; from pathlib import Path; Path(sys.argv[2]).write_text('NaN')", "error"),
])
def test_custom_execution_failure_is_not_a_success(workspace: ScientificWorkspace, code, status) -> None:
    (workspace.store.root / "failure.py").write_text(code, "utf-8")
    event = workspace.run_python("failure.py", {}, timeout_seconds=1)
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == status


def test_custom_execution_rejects_paths_outside_investigation(workspace: ScientificWorkspace) -> None:
    with pytest.raises(ValueError, match="inside this investigation"):
        workspace.run_python("../escape.py", {})


def test_custom_execution_cancellation_and_input_mutation_are_not_completed(workspace: ScientificWorkspace) -> None:
    (workspace.store.root / "mutate.py").write_text(
        "import sys; from pathlib import Path; Path(sys.argv[1]).write_text('{}'); Path(sys.argv[2]).write_text('{}')", "utf-8")
    cancelled = Event()
    cancelled.set()
    event = workspace.run_python("mutate.py", {}, cancel=cancelled)
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "cancelled"
    event = workspace.run_python("mutate.py", {})
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "error"
    assert "changed its recorded script or inputs" in record["error"]["message"]


def test_custom_execution_snapshots_caller_owned_parameters(workspace: ScientificWorkspace, monkeypatch) -> None:
    (workspace.store.root / "copy.py").write_text(
        "import sys,json; from pathlib import Path; x=json.loads(Path(sys.argv[1]).read_text()); Path(sys.argv[2]).write_text(json.dumps(x['parameters']))", "utf-8")
    parameters = {"values": [1, 2]}
    original = subprocess.Popen

    def start(*args, **kwargs):
        process = original(*args, **kwargs)
        parameters["values"].append(3)
        return process

    monkeypatch.setattr(subprocess, "Popen", start)
    event = workspace.run_python("copy.py", parameters)
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "completed"
    assert record["inputs"]["parameters"] == record["result"] == {"values": [1, 2]}
