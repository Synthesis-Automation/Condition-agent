"""Source-bound route precedent inspection, rendering, and citation regressions."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from app.web_api.scientific_presentation import present_conversation
from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.answer_finalization import _complete_empty_fields
from chem_coworker.scientific_workspace.baseline import artifact_identity, code_manifest, environment_versions
from chem_coworker.scientific_workspace.step_precedents import answer_step_precedents
from core_retrosynthesis.generic_library import build_generic_library, save_generic_library
from tests.core_retrosynthesis_tests.test_external_proposal_admission import _row, FIRST_REACTION


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def library():
    return build_generic_library(tuple(_row(FIRST_REACTION, i) for i in range(1, 5)), levels=("L0", "L1", "L2"))


@pytest.fixture
def workspace(tmp_path, library):
    paths = {"retro_library": tmp_path / "library.json.gz", "reference_catalog": tmp_path / "refs.jsonl",
             "procedure_catalog": tmp_path / "procedures.jsonl"}
    save_generic_library(library, paths["retro_library"])
    paths["reference_catalog"].write_text("\n".join(json.dumps({
        "reference_id": f"reference-{i}", "patent_number": f"US1234{i}B2", "normalized_citation": f"Patent fixture {i}",
    }) for i in range(1, 5)), "utf-8")
    paths["procedure_catalog"].write_text(json.dumps({"reaction_id": "reaction-1", "observation_id": "obs-a",
                                                   "procedure_text": "Fixture source procedure, no claimed yield."}), "utf-8")
    baseline = {"repository": str(ROOT), "code_files": code_manifest(ROOT), "environment": environment_versions(),
                "artifacts": {name: artifact_identity(path) for name, path in paths.items()}}
    InvestigationStore.create(tmp_path / ("a" * 32), objective="Inspect source precedents", baseline=baseline)
    return ScientificWorkspace(tmp_path / ("a" * 32))


def run(workspace, operation, **arguments):
    event = workspace.run(operation, arguments)
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "completed", record
    return event, record["result"]


def disconnection(workspace):
    source, result = run(workspace, "disconnect_target", target_smiles="CCN")
    return source, result["strategies"][0]["representative"]


def draft_for(record, reference):
    return _complete_empty_fields({
        "answer_markdown": "A proposal supported by template examples, not a validated preparation.",
        "molecules": [{"id": "a", "name": "Precursors", "smiles": record["selection"]["precursor_smiles"], "basis": "proposed"},
                      {"id": "b", "name": "Target", "smiles": record["selection"]["target_smiles"], "basis": "input"}],
        "steps": [{"id": "s1", "title": "Proposed step", "reactant_ids": ["a"], "product_ids": ["b"],
                   "basis": "proposed", "precedent_refs": [reference]}],
    })


def test_real_selection_page_sources_summary_and_replay(workspace):
    source, selected = disconnection(workspace)
    event, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                        realization_id=selected["realization_id"], limit=2)
    assert source.artifact_ref in event.evidence_refs
    assert record["page"]["next_offset"] == 2
    assert record["saved_match_count"] == 4
    assert record["distinct_references_on_page"] == 2
    first = record["precedents"][0]
    assert first["precursor_smiles"] and first["product_smiles"]
    assert first["reference_record"]["patent_number"] == "US12341B2"
    assert first["procedures"][0]["observation_id"] == "obs-a"
    assert first["observations"] == []
    assert first["product_comparison"]["same_constitution"] is True
    assert first["support_kind"] == "template_precedent"
    brief = workspace.call_summary(event)["result_summary"]
    assert brief["precedents"][0]["reaction_smiles"] == first["reaction_smiles"]
    assert brief["page"]["next_offset"] == 2
    _, second = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                    realization_id=selected["realization_id"], offset=2, limit=2)
    assert second["page"]["next_offset"] is None
    assert {r["reaction_id"] for r in record["precedents"]}.isdisjoint(r["reaction_id"] for r in second["precedents"])
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_assessed_step_and_route_preserve_support_and_unknown(workspace):
    source, original = run(workspace, "assess_route_proposal", proposal={"target_smiles": "CCN", "steps": [
        {"external_step_id": "amine", "target_smiles": "CCN", "precursor_smiles": "CC=O.N"},
    ]})
    _, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref, step_id="amine")
    assert record["scope"] == "saved_assessment_matches"
    assert record["precedents"]
    assert record["available_template_records"] is None
    assert workspace.store.read_artifact(source.artifact_ref)["result"] == original
    unknown, _ = run(workspace, "assess_route_step", proposal={"target_smiles": "CCN", "precursor_smiles": "CC.CN"})
    _, gap = run(workspace, "inspect_step_precedents", source_ref=unknown.artifact_ref)
    assert gap["status"] == "no_precedents_retrieved"
    assert gap["precedents"] == []
    assert gap["experimental_feasibility"] == "not_established"


@pytest.mark.parametrize("arguments", [{"realization_id": "invented"}, {"step_id": "invented"},
                                       {"offset": -1}, {"limit": 100}])
def test_invalid_selection_is_recorded_as_error(workspace, arguments):
    source, _ = disconnection(workspace)
    event = workspace.run("inspect_step_precedents", {"source_ref": source.artifact_ref, **arguments})
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "error"


def test_answer_binding_rendering_and_tamper_recovery(workspace):
    source, selected = disconnection(workspace)
    event, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                        realization_id=selected["realization_id"], limit=1)
    draft = draft_for(record, event.artifact_ref)
    answer = ScientificAnswer.model_validate(draft)
    assert event.artifact_ref in validate_answer_evidence(answer, workspace.store)
    evidence = answer_step_precedents(workspace.store, draft)
    conversation = {"id": "a" * 32, "turns": [{"question": "Prepare target", "answer": draft,
                                                "step_precedent_evidence": evidence}]}
    before = deepcopy(conversation)
    view = present_conversation(conversation)["turns"][0]["structured_presentation"]
    first = view["steps"][0]["supporting_evidence"][0]["precedents"][0]
    assert first["drawing_status"] == "drawn"
    assert first["reference_url"] == "https://patents.google.com/patent/US12341B2"
    assert conversation == before
    changed = deepcopy(draft)
    changed["molecules"][1]["smiles"] = "CCO"
    with pytest.raises(ValueError, match="does not match"):
        validate_answer_evidence(ScientificAnswer.model_validate(changed), workspace.store)
    changed["molecules"][1]["smiles"] = "NCC"
    validate_answer_evidence(ScientificAnswer.model_validate(changed), workspace.store)
    artifact = workspace.store.root / "artifacts" / (event.artifact_ref.split(":")[1] + ".json")
    artifact.write_text("{}", "utf-8")
    assert answer_step_precedents(workspace.store, draft)["s1"][0]["status"] == "evidence_unavailable"


def test_condition_observations_remain_separate(workspace, monkeypatch):
    source, selected = disconnection(workspace)
    monkeypatch.setattr(workspace.operations, "get_precedents", lambda *a, **kw: {"records": [
        {"reaction_id": "reaction-1", "observation_id": "obs-a", "yield_pct": 40},
        {"reaction_id": "reaction-1", "observation_id": "obs-b", "yield_pct": 80},
        {"reaction_id": "unrelated", "observation_id": "bad", "yield_pct": 99},
    ]})
    _, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                    realization_id=selected["realization_id"], limit=1)
    first = record["precedents"][0]
    assert [row["yield_pct"] for row in first["observations"]] == [40, 80]
    assert "yield_pct" not in first
    assert first["experimental_link_scope"] == "reaction_id_only_template_has_no_observation_id"


def test_unrelated_call_and_wrong_stereochemistry_cannot_support_step(workspace):
    source, _ = run(workspace, "analyze_molecule", smiles="CCO")
    invalid = workspace.run("inspect_step_precedents", {"source_ref": source.artifact_ref})
    assert workspace.store.read_artifact(invalid.artifact_ref)["execution_status"] == "error"
    from chem_coworker.scientific_workspace.step_precedents import load_step_precedent_evidence, SCHEMA_VERSION

    record = workspace.store.append("call", {"operation": "inspect_step_precedents", "execution_status": "completed",
        "result": {"schema_version": SCHEMA_VERSION, "selection": {
            "precursor_smiles": "CCO", "target_smiles": "N[C@H](C)C(=O)O"}}})
    with pytest.raises(ValueError, match="does not match"):
        load_step_precedent_evidence(workspace.store, record.artifact_ref, "CCO", "N[C@@H](C)C(=O)O")


def test_missing_reference_catalog_remains_a_gap(workspace, monkeypatch):
    source, selected = disconnection(workspace)
    original = workspace.operations._path

    def path(name):
        if name == "reference_catalog":
            raise FileNotFoundError("Not configured")
        return original(name)

    monkeypatch.setattr(workspace.operations, "_path", path)
    _, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                    realization_id=selected["realization_id"], limit=1)
    assert record["precedents"] and record["precedents"][0]["reference_record"] is None
    assert record["reference_catalog_status"] == "catalog_unavailable"


def test_saved_conversation_api_resolves_precedents_without_running_science(workspace, monkeypatch):
    from fastapi.testclient import TestClient
    from app.web_api.main import create_app
    from chem_coworker.scientific_workspace.conversation import ConversationService

    source, selected = disconnection(workspace)
    event, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                        realization_id=selected["realization_id"], limit=1)
    draft = draft_for(record, event.artifact_ref)
    identity, turn_id = "a" * 32, "b" * 32
    (workspace.store.root / "conversation.json").write_text(json.dumps({"id": identity, "title": "Review", "created_at": "2026-09-29"}), "utf-8")
    directory = workspace.store.root / "turns" / turn_id
    directory.mkdir(parents=True)
    saved = {"id": turn_id, "created_at": "2026-09-29", "status": "completed", "question": "Prepare this", "answer": draft,
             "activity_version": 4, "progress": []}
    path = directory / "turn.json"
    path.write_text(json.dumps(saved), "utf-8")
    original = path.read_bytes()
    monkeypatch.setattr(ScientificWorkspace, "run", lambda *a, **kw: pytest.fail("Read path ran science"))
    service = ConversationService(workspace.store.root.parent, runtime=object())
    try:
        client = TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False), base_url="http://127.0.0.1")
        response = client.get(f"/api/v1/scientific/conversations/{identity}")
        assert response.status_code == 200
        evidence = response.json()["turns"][0]["structured_presentation"]["steps"][0]["supporting_evidence"]
        assert evidence[0]["precedents"][0]["drawing_status"] == "drawn"
        assert path.read_bytes() == original
    finally:
        service.close()


def test_final_assessment_requires_inspection_and_preserves_old_answer_preview(workspace, monkeypatch):
    source, assessed = run(workspace, "assess_route_proposal", proposal={"target_smiles": "CCN", "steps": [
        {"external_step_id": "chosen", "target_smiles": "CCN", "precursor_smiles": "CC=O.N"},
    ]})
    draft = draft_for({"selection": assessed["proposal"]["steps"][0]}, source.artifact_ref)
    draft["steps"][0]["precedent_refs"] = []
    draft["evidence_refs"] = [source.artifact_ref]
    before = deepcopy(draft)
    path = workspace.store.root / "turns" / "test" / "answer-draft.json"
    path.parent.mkdir(parents=True)
    with pytest.raises(ValueError, match='inspect_step_precedents.*chosen'):
        workspace.finalize_answer(path, draft)
    assert not path.exists()
    with pytest.raises(ValueError, match="Available supporting reactions"):
        validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)
    # A read of old answers must not rerun science or silently claim agent review.
    with monkeypatch.context() as patch:
        patch.setattr(ScientificWorkspace, "run", lambda *a, **kw: pytest.fail("Read ran science"))
        view = answer_step_precedents(workspace.store, draft)["s1"][0]
    assert view["artifact_ref"] == source.artifact_ref
    assert view["evidence_origin"] == "saved_assessment"
    assert view["inspection_status"] == "not_attached"
    assert view["precedents"][0]["reaction_smiles"] == "CC=O.N>>CCN"
    assert view["precedents"][0]["observations"] == []
    assert draft == before
    inspection, _ = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref, step_id="chosen")
    draft["steps"][0]["precedent_refs"] = [inspection.artifact_ref]
    workspace.finalize_answer(path, draft)
    assert json.loads(path.read_text("utf-8"))["steps"][0]["precedent_refs"] == [inspection.artifact_ref]


def test_available_inspection_must_be_linked_and_can_be_recovered_for_display(workspace):
    source, selected = disconnection(workspace)
    event, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                        realization_id=selected["realization_id"])
    draft = draft_for(record, event.artifact_ref)
    draft["steps"][0]["precedent_refs"] = []
    with pytest.raises(ValueError, match="attach existing inspection"):
        validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)
    recovered = answer_step_precedents(workspace.store, draft)["s1"][0]
    assert recovered["artifact_ref"] == event.artifact_ref
    assert recovered["evidence_origin"] == "recovered_inspection"
    assert recovered["precedents"] == record["precedents"]
    # An earlier route to the same product with different reactants is not support.
    draft["molecules"][0]["smiles"] = "CC.CN"
    validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)
    assert answer_step_precedents(workspace.store, draft)["s1"] == []


def test_empty_inspection_cannot_hide_known_matches(workspace):
    source, selected = disconnection(workspace)
    event, record = run(workspace, "inspect_step_precedents", source_ref=source.artifact_ref,
                        realization_id=selected["realization_id"])
    empty = deepcopy(record)
    empty.update(precedents=[], saved_match_count=0, status="no_precedents_retrieved")
    saved = workspace.store.append("call", {"operation": "inspect_step_precedents", "execution_status": "completed", "result": empty})
    draft = draft_for(record, saved.artifact_ref)
    with pytest.raises(ValueError, match="Available supporting reactions"):
        validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)


def test_recovery_preserves_stereo_and_does_not_use_corrupt_sources(workspace):
    source, assessed = run(workspace, "assess_route_step", proposal={"target_smiles": "CCN", "precursor_smiles": "CC=O.N"})
    # Exercise binding without requiring a new chemistry operator fixture.
    payload = {"operation": "assess_route_step", "execution_status": "completed", "result": deepcopy(assessed)}
    payload["result"]["proposal"]["target_smiles"] = "N[C@H](C)C(=O)O"
    saved = workspace.store.append("call", payload)
    draft = draft_for({"selection": payload["result"]["proposal"]}, saved.artifact_ref)
    draft["steps"][0]["precedent_refs"] = []
    assert answer_step_precedents(workspace.store, draft)["s1"]
    draft["molecules"][1]["smiles"] = "N[C@@H](C)C(=O)O"
    assert answer_step_precedents(workspace.store, draft)["s1"] == []
    validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)
    draft["molecules"][1]["smiles"] = "N[C@H](C)C(=O)O"
    (workspace.store.root / "artifacts" / (saved.artifact_ref.split(":")[1] + ".json")).write_text("{}", "utf-8")
    assert answer_step_precedents(workspace.store, draft)["s1"][0]["status"] == "evidence_unavailable"
