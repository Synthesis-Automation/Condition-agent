"""Regressions for exact candidate selection, source leads and recipe association."""

from copy import deepcopy
from dataclasses import asdict

import pytest

from chem_coworker.scientific_workspace.adapters.step_selection import _selection
from chem_coworker.scientific_workspace.answers.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.answers.answer_finalization import _complete_empty_fields
from chem_coworker.scientific_workspace.views.step_assessments import answer_step_assessments
from condition_recommender import assess_reaction_recipe
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from tests.chem_coworker.test_route_investigation_workspace import library, workspace, call


def ambiguous_disconnection():
    """Two actual mixed-ester Hantzsch inputs share a template realization ID."""
    target = "CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccc(F)c(F)c1"
    inputs = ["CCOC(=O)CC(C)=O.COC(=O)/C=C(/C)N.O=Cc1ccc(F)c(F)c1",
              "CCOC(=O)/C=C(/C)N.COC(=O)CC(C)=O.O=Cc1ccc(F)c(F)c1"]
    return {"operation": "disconnect_target", "execution_status": "completed", "result": {"strategies": [
        {"strategy_id": f"s{i}", "representative": {"realization_id": "REAL2:shared", "target_smiles": target,
         "precursor_smiles": precursor, "operator_id": "OP1:same", "template_id": "GRT3:same"}}
        for i, precursor in enumerate(inputs)]}}


def test_shared_template_id_never_silently_selects_different_precursors():
    payload = ambiguous_disconnection()
    before = deepcopy(payload)
    with pytest.raises(ValueError, match="Ambiguous"):
        _selection(payload, None, "REAL2:shared")
    selected, _ = _selection(payload, None, "REAL2:shared", strategy_id="s1")
    assert selected["precursor_smiles"] == payload["result"]["strategies"][1]["representative"]["precursor_smiles"]
    reordered = ".".join(reversed(selected["precursor_smiles"].split(".")))
    assert _selection(payload, None, "REAL2:shared", precursor_smiles=reordered)[0] == selected
    with pytest.raises(ValueError, match="No saved candidate"):
        _selection(payload, None, "REAL2:shared", strategy_id="s0", precursor_smiles=reordered)
    assert payload == before


def test_same_graph_duplicates_are_allowed_but_stereoisomers_remain_distinct():
    payload = ambiguous_disconnection()
    first, second = [x["representative"] for x in payload["result"]["strategies"]]
    second["precursor_smiles"] = ".".join(reversed(first["precursor_smiles"].split(".")))
    assert _selection(payload, None, "REAL2:shared")[0] == first
    first["precursor_smiles"] = "N[C@H](C)C(=O)O"
    second["precursor_smiles"] = "N[C@@H](C)C(=O)O"
    with pytest.raises(ValueError, match="Ambiguous"):
        _selection(payload, None, "REAL2:shared")
    assert _selection(payload, None, "REAL2:shared", precursor_smiles=second["precursor_smiles"])[0] == second


def test_ambiguous_selection_errors_are_recorded_for_inspection_and_validity(workspace):
    source = workspace.store.append("call", ambiguous_disconnection())
    for operation in ("inspect_step_precedents", "assess_retro_validity"):
        event = workspace.run(operation, {"source_ref": source.artifact_ref, "realization_id": "REAL2:shared"})
        saved = workspace.store.read_artifact(event.artifact_ref)
        assert saved["execution_status"] == "error"
        assert "Ambiguous" in saved["error"]["message"]


def test_upstream_search_checks_all_roots_and_retains_exact_excerpt_offsets(workspace):
    source = workspace.capture_source("Upstream: acetaldehyde was prepared from ethanol.\nLater acetaldehyde was used.",
                                      url="https://example.org/patent")
    excerpt = workspace.record_source_excerpt(source.artifact_ref, excerpt="Later acetaldehyde was used.")
    failed = workspace.store.append("literature_source", {"extraction": {"text": ""}, "retrieval_status": "failed"})
    event, result = call(workspace, "search_captured_sources", queries=["acetaldehyde"],
                         source_refs=[excerpt.artifact_ref, source.artifact_ref], limit=1)
    assert result["total"] == 2 and result["next_offset"] == 1
    assert result["matches"][0]["match_start"] < 20
    match = result["matches"][0]
    recorded = workspace.record_source_excerpt(match["source_ref"], start=match["start"], end=match["end"])
    assert workspace.store.read_artifact(recorded.artifact_ref)["text"] == match["text"]
    assert source.artifact_ref in event.evidence_refs
    _, default = call(workspace, "search_captured_sources", queries=["missing"])
    assert default["total"] == 0 and failed.artifact_ref not in default["evidence_refs"]


def test_long_source_search_pages_interleaved_terms_without_losing_counts(workspace):
    source = workspace.capture_source("Compound 1. " * 20000, url="https://example.org/long-patent")
    _, result = call(workspace, "search_captured_sources", queries=["Compound", "1", "Compound"],
                     source_refs=[source.artifact_ref], offset=39998, limit=2)
    assert result["total"] == 40000 and result["next_offset"] is None
    assert [(item["query"], item["match_start"]) for item in result["matches"]] == [
        ("Compound", 19999 * 12), ("1", 19999 * 12 + 9)]
    assert all(len(item["text"]) <= 607 for item in result["matches"])


def test_route_leaf_report_retains_assumed_terminal_and_unsearched_leaves(workspace):
    source = workspace.capture_source("Preparation of acetaldehyde: oxidation of ethanol.", url="https://example.org/patent")
    route, _ = call(workspace, "assess_route_proposal", proposal={"target_smiles": "CCN", "steps": [
        {"external_step_id": "s1", "target_smiles": "CCN", "precursor_smiles": "CC=O.N"}]})
    event, result = call(workspace, "inspect_route_inputs", source_ref=route.artifact_ref,
                         leaf_queries=[{"smiles": "CC=O", "terms": ["acetaldehyde"]}])
    aldehyde = next(item for item in result["leaves"] if item["smiles"] == "CC=O")
    assert aldehyde["captured_source_search"]["total"] == 1
    assert aldehyde["starting_material_assessment"]["status"] == "assumed_terminal"
    assert any(item["source_search_status"] == "terms_not_supplied" for item in result["leaves"])
    assert source.artifact_ref in event.evidence_refs
    invalid = workspace.run("inspect_route_inputs", {"source_ref": route.artifact_ref,
                            "leaf_queries": [{"smiles": "CCO", "terms": ["ethanol"]}]})
    assert workspace.store.read_artifact(invalid.artifact_ref)["execution_status"] == "error"


def proposed_draft():
    return _complete_empty_fields({"answer_markdown": "Proposed reaction; evidence is incomplete.",
        "molecules": [{"id": "r", "name": "Inputs", "smiles": "CCBr.N", "basis": "proposed"},
                      {"id": "p", "name": "Product", "smiles": "CCN", "basis": "proposed"}],
        "steps": [{"id": "s1", "title": "Substitution", "basis": "proposed", "reactant_ids": ["r"],
                   "product_ids": ["p"], "conditions": [{"text": "Proposed ethanol medium", "basis": "proposed"}]}]})


def test_actual_recipe_is_normalized_and_bound_to_answer_not_an_analogue(workspace):
    event, result = call(workspace, "assess_proposed_recipe", reaction_smiles="CCBr.N>>CCN",
        components=[{"raw_identifier": "ethanol", "source_field": "proposed", "identifier_type": "name",
                     "source_role_hint": "solvent"}],
        operating_conditions={"stages": [{"stage_index": 0, "temperature_c": 0, "time_h": 0.25},
                                         {"stage_index": 1, "temperature_c": 25, "time_h": 2}]})
    assert result["proposed_recipe"]["solvents"][0]["identity_status"] == "resolved"
    assert len(result["proposed_recipe"]["stages"]) == 2
    assert canonical_bytes(result["compatibility"]) == canonical_bytes(asdict(assess_reaction_recipe("CCBr.N>>CCN", result["proposed_recipe"])))
    draft = proposed_draft()
    assert any(x["gap"] == "actual_recipe_check_not_attached" for x in workspace.answer_preflight(draft)["warnings"])
    attached = workspace.attach_recipe_check(draft, "s1", event.artifact_ref)
    assert "recipe_assessment_refs" not in draft["steps"][0]
    assert event.artifact_ref in validate_answer_evidence(ScientificAnswer.model_validate(attached), workspace.store)
    view = answer_step_assessments(workspace.store, attached)["s1"]
    assert view["recipe_assessments"][0]["status"] == result["compatibility"]["status"]
    attached["molecules"][0]["smiles"] = "CCCl.N"
    with pytest.raises(ValueError, match="does not match"):
        validate_answer_evidence(ScientificAnswer.model_validate(attached), workspace.store)


def test_unresolved_recipe_check_cannot_claim_conflict_or_feasibility(workspace):
    _, result = call(workspace, "assess_proposed_recipe", reaction_smiles="CC.CN>>CCN",
                     components=[{"raw_identifier": "unknown fixture substance", "source_field": "proposal"}],
                     operating_conditions={})
    assert result["compatibility"]["status"] == "unknown"
    assert result["experimental_feasibility"] == "not_established"
    assert not result["compatibility"]["hard_conflicts"]


def test_preflight_keeps_leaf_assumptions_visible_after_input_inspection(workspace):
    route, _ = call(workspace, "assess_route_proposal", proposal={"target_smiles": "CCN", "steps": [
        {"external_step_id": "s1", "target_smiles": "CCN", "precursor_smiles": "CCBr.N"}]})
    inspection, _ = call(workspace, "inspect_step_precedents", source_ref=route.artifact_ref, step_id="s1")
    draft = proposed_draft()
    draft["steps"][0]["precedent_refs"] = [inspection.artifact_ref]
    # Preflight follows this inspection's parent route even without a direct citation.
    first = workspace.answer_preflight(draft)
    assert first["valid"] and any(item["gap"] == "route_inputs_not_inspected" for item in first["warnings"])
    inputs, _ = call(workspace, "inspect_route_inputs", source_ref=route.artifact_ref, leaf_queries=[])
    second = workspace.answer_preflight(draft)
    assert second["valid"] and any(item.get("artifact_ref") == inputs.artifact_ref for item in second["warnings"])
    assert not any(item["gap"] == "route_inputs_not_inspected" for item in second["warnings"])


@pytest.mark.parametrize("conditions", [{"unsupported": 1}, {"temperature_c": "90"},
    {"stages": [{"stage_index": 0}, {"stage_index": 0}]}])
def test_recipe_operating_input_does_not_silently_drop_invalid_fields(workspace, conditions):
    event = workspace.run("assess_proposed_recipe", {"reaction_smiles": "CCBr.N>>CCN",
        "components": [{"raw_identifier": "water", "source_field": "proposal"}], "operating_conditions": conditions})
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "error"
