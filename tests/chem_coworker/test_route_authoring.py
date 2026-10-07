"""Compact route authoring preserves topology, provenance and publication gates."""

from copy import deepcopy
import base64
import json
import xml.etree.ElementTree as ET

import pytest

from app.web_api.scientific_presentation import present_conversation
from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.adapters.step_selection import SCHEMA_VERSION
from chem_coworker.scientific_workspace.answers.answer_contracts import ScientificAnswer
from chem_coworker.scientific_workspace.core.store import InvestigationStore


@pytest.fixture
def workspace(tmp_path):
    root = tmp_path / "investigation"
    InvestigationStore.create(root, objective="Compact route transport", baseline={})
    attempt = root / "turns" / "attempt"
    attempt.mkdir(parents=True)
    return ScientificWorkspace(root), attempt / "answer-draft.json"


def proposal():
    return {"target_smiles": "CC(=O)O", "routes": [{"steps": [
        {"reaction_smiles": "OCC>>CC=O"},
        {"reaction_smiles": "O=CC>>CC(=O)O"},
    ]}]}


def test_minimal_route_finalizes_without_review_or_scientific_calls_and_draws_svg(workspace):
    w, path = workspace
    request = proposal()
    original = deepcopy(request)
    draft = w.route_answer(request)
    receipt = w.finalize_answer(path, draft)
    assert receipt["answer_file"] == "answer-draft.json"
    assert request == original
    assert not w.store.events()
    saved = json.loads(path.read_text("utf-8"))
    answer = ScientificAnswer.model_validate(saved)
    assert len(answer.molecules) == 3
    assert answer.steps[1].after_step_ids == [answer.steps[0].id]
    assert answer.steps[1].reactant_ids == answer.steps[0].product_ids
    assert all(step.basis == "proposed" and not step.conditions for step in answer.steps)
    assert all(step.rationale is None and not step.literature_reactions for step in answer.steps)
    assert not answer.claims and not answer.sources
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": saved}]})
    structured = view["turns"][0]["structured_presentation"]
    assert structured["routes"][0]["drawing_status"] == "linear_scheme"
    assert structured["steps"][0]["reaction_smiles"] == "CCO>>CC=O"
    svg = ET.fromstring(base64.b64decode(structured["routes"][0]["image_url"].split(",", 1)[1]))
    assert svg.tag.endswith("svg")


def test_literature_metadata_and_cautions_are_reused_without_source_drawings(workspace):
    w, path = workspace
    source = w.capture_source("Example 7: reported on an analogue, not the proposed target.",
                              url="https://example.org/paper", title="Paper", locator="Example 7")
    support = {"ref": source.artifact_ref, "locator": "Example 7", "note": "Different substrate."}
    request = proposal()
    for step in request["routes"][0]["steps"]:
        step.update(conditions="Proposed screening conditions", support=[support])
    draft = w.route_answer(request)
    w.finalize_answer(path, draft)
    saved = json.loads(path.read_text("utf-8"))
    assert len(saved["sources"]) == 1
    assert saved["sources"][0]["url"] == "https://example.org/paper"
    assert saved["sources"][0]["title"] == "Paper"
    assert saved["evidence_refs"] == [source.artifact_ref]
    assert saved["steps"][0]["conditions"][0]["basis"] == "proposed"
    assert saved["steps"][0]["limitations"] == ["Example 7: Different substrate."]
    assert not saved["steps"][0]["literature_reactions"]
    assert [event.kind for event in w.store.events()] == ["literature_source"]


def test_route_local_branches_and_explicit_ambiguous_dependencies(workspace):
    w, _ = workspace
    request = {"target_smiles": "CCN", "routes": [{"steps": [
        {"reaction_smiles": "CCO>>CC=O"}, {"reaction_smiles": "[NH4+]>>N"},
        {"reaction_smiles": "CC=O.N>>CCN"},
    ]}, {"steps": [{"reaction_smiles": "CC=O.N>>CCN"}]}]}
    draft = w.route_answer(request)
    assert draft["steps"][2]["after_step_ids"] == ["r1_s1", "r1_s2"]
    assert draft["steps"][3]["after_step_ids"] == []
    request["routes"][0]["steps"].insert(1, {"reaction_smiles": "CC(O)O>>CC=O"})
    with pytest.raises(ValueError, match="Ambiguous intermediate"):
        w.route_answer(request)
    request["routes"][0]["steps"][-1]["after_steps"] = [2, 3]
    assert w.route_answer(request)["steps"][3]["after_step_ids"] == ["r1_s2", "r1_s3"]


def test_stereo_and_component_multiplicity_are_preserved(workspace):
    w, _ = workspace
    request = {"target_smiles": "CCN", "routes": [{"steps": [
        {"reaction_smiles": "CCO>>C[C@H](O)Cl"},
        {"reaction_smiles": "C[C@@H](O)Cl.N.N>>CCN"},
    ]}]}
    draft = w.route_answer(request)
    assert draft["steps"][1]["after_step_ids"] == []
    assert len(draft["steps"][1]["reactant_ids"]) == 3
    assert draft["steps"][1]["reactant_ids"][-1] == draft["steps"][1]["reactant_ids"][-2]


@pytest.mark.parametrize("step", [
    {"reaction_smiles": "CCO>N>CC=O"},
    {"reaction_smiles": "CCO>>"},
    {"reaction_smiles": "C1>>CC=O"},
    {"reaction_smiles": "CCO..N>>CC=O"},
    {"reaction_smiles": "CCO>>CC=O", "after_steps": [1]},
    {"reaction_smiles": "CCO>>CC=O", "basis": "reported"},
    {"reaction_smiles": "CCO>>CC=O", "conditions": "Reported", "conditions_basis": "reported"},
    {"reaction_smiles": "CCO>>CC=O", "invented_field": "reject"},
])
def test_invalid_input_or_unsupported_reported_claim_is_not_repaired(workspace, step):
    w, path = workspace
    with pytest.raises(ValueError):
        w.route_answer({"target_smiles": "CC=O", "routes": [{"steps": [step]}]})
    assert not path.exists() and not w.store.events()


def test_incomplete_route_remains_incomplete(workspace):
    w, _ = workspace
    request = proposal()
    request["routes"][0]["steps"].pop()
    draft = w.route_answer(request)
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": draft}]})
    assert view["turns"][0]["structured_presentation"]["routes"][0]["unreached_target_ids"]


def test_saved_inspections_attach_and_exact_step_publication_rules_still_apply(workspace):
    w, path = workspace
    inspection = w.store.append("call", {
        "operation": "inspect_step_precedents", "execution_status": "completed",
        "result": {"schema_version": SCHEMA_VERSION,
                   "selection": {"precursor_smiles": "CCO", "target_smiles": "CC=O"},
                   "precedents": [{"reaction_id": "fixture"}]},
    })
    request = {"target_smiles": "CC=O", "routes": [{"steps": [{"reaction_smiles": "OCC>>O=CC"}]}]}
    with pytest.raises(ValueError, match="Available supporting reactions"):
        w.finalize_answer(path, w.route_answer(request))
    request["routes"][0]["steps"][0]["support"] = [{"ref": inspection.artifact_ref, "locator": "result.precedents"}]
    draft = w.route_answer(request)
    assert draft["steps"][0]["precedent_refs"] == [inspection.artifact_ref]
    w.finalize_answer(path, draft)
    before = path.read_bytes()
    request["routes"][0]["steps"][0]["reaction_smiles"] = "CCO>>CC(=O)O"
    with pytest.raises(ValueError, match="does not match"):
        w.finalize_answer(path, w.route_answer(request))
    assert path.read_bytes() == before


def test_failed_scientific_source_cannot_become_route_support(workspace):
    w, _ = workspace
    failure = w.store.append("call", {"operation": "analyze_molecule", "execution_status": "error"})
    request = proposal()
    request["routes"][0]["steps"][0]["support"] = [{"ref": failure.artifact_ref, "locator": "result"}]
    with pytest.raises(ValueError, match="completed recorded"):
        w.route_answer(request)


def test_route_input_help_is_available_without_loading_full_answer_schema(workspace):
    w, _ = workspace
    help_entries = w.help(["route_answer", "inspect_artifact"])
    assert help_entries[0]["input_schema"]["properties"]["schema_version"]["const"] == "route_answer_input.v1"
    assert "inspect_artifact" in help_entries[1]["example"]
