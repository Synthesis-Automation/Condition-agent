"""Saved assessment attribution: exact sides, uncertainty, history and integrity."""

from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.web_api.main import create_app
from app.web_api.scientific_presentation import present_conversation
from chem_coworker.scientific_workspace.core.store import InvestigationStore
from chem_coworker.scientific_workspace.runtime.activity import ACTIVITY_VERSION
from chem_coworker.scientific_workspace.runtime.conversation import ConversationService
from chem_coworker.scientific_workspace.views.step_assessments import answer_step_assessments
from condition_recommender import assess_reaction_recipe


@pytest.fixture
def saved(tmp_path: Path) -> tuple[InvestigationStore, dict]:
    store = InvestigationStore.create(tmp_path / ("a" * 32), objective="View regression", baseline={})
    proposed = {"basis": "proposed", "source_ids": [], "limitations": []}
    answer = {
        "schema_version": "scientific_answer.v2", "answer_markdown": "Proposed reaction.",
        "evidence_refs": [], "uncertainties": [], "needs_user_input": False, "sources": [],
        "molecules": [{"id": "r", "name": "Inputs", "smiles": "CCBr.N", **proposed},
                      {"id": "p", "name": "Product", "smiles": "CCN", **proposed}],
        "target_molecule_ids": ["p"], "steps": [{
            "id": "s1", "title": "Proposed substitution", "reactant_ids": ["r"], "product_ids": ["p"],
            "after_step_ids": [], "conditions": [], "yield_info": None, **proposed,
        }], "routes": [], "claims": [],
    }
    return store, answer


def _structural(store: InvestigationStore, precursors: str, status: str = "precedent_supported") -> str:
    """Save an assessment-shaped attribution fixture, not a computed chemistry claim."""
    return store.append("call", {
        "operation": "assess_route_step", "execution_status": "completed",
        "result": {"proposal": {"precursor_smiles": precursors, "target_smiles": "CCN"},
                   "assessment": {"status": status, "actionable": False, "admission_eligible": False,
                                  "gates": [{"gate_id": "condition_support", "status": "not_run",
                                             "summary": "No condition evaluation requested.", "warnings": []}]}},
    }).artifact_ref


def test_exact_structure_matching_does_not_hide_conflicting_assessments(saved) -> None:
    store, answer = saved
    unrelated = _structural(store, "CCCl.N")
    uncited = _structural(store, "CCBr.N", "conflicting")
    first = _structural(store, "N.CCBr")
    second = _structural(store, "CCBr.N", "ambiguous")
    answer["evidence_refs"] = [unrelated, first, second]
    view = answer_step_assessments(store, answer)["s1"]
    assert view["status"] == "recorded"
    assert {r["artifact_ref"] for r in view["structural_assessments"]} == {first, second, uncited}
    assert next(r for r in view["structural_assessments"] if r["artifact_ref"] == uncited)["attribution"] == "saved_exact_step_check"
    assert {r["status"] for r in view["structural_assessments"]} == {"precedent_supported", "ambiguous", "conflicting"}
    assert not view["recipe_assessments"]
    assert all(r["gates"][0]["status"] == "not_run" for r in view["structural_assessments"])


def test_uncited_route_failure_is_visible_but_later_checks_and_recipes_are_excluded(saved):
    store, answer = saved
    before = deepcopy(answer)
    proposal = {"external_step_id": "original", "precursor_smiles": "N.CCBr", "target_smiles": "CCN"}
    failed = store.append("call", {
        "operation": "assess_route_proposal", "execution_status": "completed",
        "result": {"proposal": {"steps": [proposal]}, "assessment": {"step_assessments": [{
            "external_step_id": "original", "assessment": {"status": "ambiguous", "admission_eligible": False,
            "gates": [{"gate_id": "atom_correspondence", "status": "unresolved"}]}}]},
            "proposed_recipe_assessments": {"original": {"status": "conflicting"}}},
    }).artifact_ref
    receipt = store.append("assistant_answer", answer)
    later = _structural(store, "CCBr.N")
    view = answer_step_assessments(store, answer, answer_ref=receipt.artifact_ref)["s1"]
    assert view["status"] == "recorded"
    assert [r["artifact_ref"] for r in view["structural_assessments"]] == [failed]
    assert view["structural_assessments"][0]["admission_eligible"] is False
    assert not view["recipe_assessments"]
    assert later not in str(view)
    assert answer == before


def test_recipe_coverage_reaches_presentation_without_modifying_saved_answer(saved) -> None:
    store, answer = saved
    recipe = {"recipe_id": "test-recipe", "solvents": [{"identity_status": "resolved"}]}
    event = store.append("call", {
        "operation": "assess_recipe", "execution_status": "completed",
        "arguments": {"reaction_smiles": "N.CCBr>>CCN", "recipe": recipe},
        "result": asdict(assess_reaction_recipe("N.CCBr>>CCN", recipe)),
    })
    answer["evidence_refs"] = [event.artifact_ref]
    original = deepcopy(answer)
    evidence = answer_step_assessments(store, answer)
    conversation = {"id": "a" * 32, "turns": [{"question": "Plan", "answer": answer,
                                                "step_assessment_evidence": evidence}]}
    view = present_conversation(conversation)["turns"][0]["structured_presentation"]
    check = view["steps"][0]["assessment_evidence"]["recipe_assessments"][0]
    assert check["coverage"]["capability_status"] == "not_covered"
    assert check["recipe_id"] == "test-recipe"
    assert check["artifact_url"].endswith(event.artifact_ref)
    assert answer == original
    assert "artifact_url" not in evidence["s1"]["recipe_assessments"][0]


def test_historical_recipe_has_unreported_coverage_and_corruption_fails_closed(saved) -> None:
    store, answer = saved
    event = store.append("call", {
        "operation": "assess_recipe", "execution_status": "completed",
        "arguments": {"reaction_smiles": "CCBr.N>>CCN", "recipe": {}},
        "result": {"status": "no_known_conflict", "compatible": True, "score": 1.0},
    })
    answer["evidence_refs"] = [event.artifact_ref]
    view = answer_step_assessments(store, answer)["s1"]
    assert view["recipe_assessments"][0]["coverage"] is None
    (store.root / "artifacts" / (event.artifact_ref.split(":")[1] + ".json")).write_text("{}")
    broken = answer_step_assessments(store, answer)["s1"]
    assert broken["status"] == "evidence_unavailable"
    assert not broken["recipe_assessments"]


def test_missing_reactant_or_different_stereo_does_not_recover_an_assessment(saved) -> None:
    store, answer = saved
    reference = _structural(store, "CCBr")
    answer["evidence_refs"] = [reference]
    assert answer_step_assessments(store, answer)["s1"]["status"] == "not_recorded"
    answer["molecules"][0]["smiles"] = "C[C@H](N)CC"
    answer["molecules"][1]["smiles"] = "C[C@H](O)CC"
    event = store.append("call", {
        "operation": "assess_recipe", "execution_status": "completed",
        "arguments": {"reaction_smiles": "C[C@@H](N)CC>>C[C@@H](O)CC", "recipe": {}},
        "result": {"status": "unknown"},
    })
    answer["evidence_refs"] = [event.artifact_ref]
    assert answer_step_assessments(store, answer)["s1"]["status"] == "not_recorded"


def test_conversation_api_exposes_saved_assessments_without_running_agent(saved) -> None:
    store, answer = saved
    answer["evidence_refs"] = [_structural(store, "N.CCBr")]
    identity, turn_id = store.root.name, "b" * 32
    (store.root / "conversation.json").write_text(json.dumps({"id": identity}), encoding="utf-8")
    directory = store.root / "turns" / turn_id
    directory.mkdir(parents=True)
    state = {"id": turn_id, "created_at": "2026-10-04T00:00:00+00:00", "question": "Plan",
             "status": "completed", "activity_version": ACTIVITY_VERSION, "progress": [], "answer": answer}
    original = json.dumps(state)
    (directory / "turn.json").write_text(original, encoding="utf-8")
    service = ConversationService(store.root.parent, runtime=object())
    try:
        client = TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False),
                            base_url="http://127.0.0.1")
        response = client.get(f"/api/v1/scientific/conversations/{identity}")
        assert response.status_code == 200
        turn = response.json()["turns"][0]
        evidence = turn["structured_presentation"]["steps"][0]["assessment_evidence"]
        assert evidence["structural_assessments"][0]["status"] == "precedent_supported"
        assert turn["step_assessment_evidence"]["s1"]["status"] == "recorded"
        assert (directory / "turn.json").read_text(encoding="utf-8") == original
    finally:
        service.close()
