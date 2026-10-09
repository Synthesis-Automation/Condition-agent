"""Saved operator mappings remain proposed evidence and pass canonical validation."""

from copy import deepcopy

import pytest

from tests.chem_coworker.test_route_investigation_workspace import library, workspace, call
from tests.chem_coworker.test_investigation_review_improvements import ambiguous_disconnection


def test_selected_mapping_survives_step_route_validity_and_revision(workspace):
    source, result = call(workspace, "disconnect_target", target_smiles="CCN", top_k=2,
                          max_templates_to_apply=10, max_candidates_to_validate=5)
    strategy = result["strategies"][0]
    candidate = strategy["representative"]
    proposal = {"target_smiles": candidate["target_smiles"], "precursor_smiles": candidate["precursor_smiles"],
                "saved_candidate": {"source_ref": source.artifact_ref,
                                    "realization_id": candidate["realization_id"],
                                    "strategy_id": strategy["strategy_id"]}}
    original = deepcopy(proposal)
    step, assessed = call(workspace, "assess_route_step", proposal=proposal)
    assert assessed["proposal"]["mapped_reaction_smiles"] == candidate["condition_query_reaction_smiles"]
    assert assessed["mapping_evidence"]["origin"] == "saved_operator_reconstruction"
    assert assessed["mapping_evidence"]["observed_reaction"] is False
    assert source.artifact_ref in step.evidence_refs
    assert proposal == original
    route, reviewed = call(workspace, "assess_route_proposal", proposal={"target_smiles": "CCN", "steps": [
        {"external_step_id": "s1", **proposal}]})
    assert reviewed["mapping_evidence"]["s1"] == assessed["mapping_evidence"]
    assert source.artifact_ref in route.evidence_refs
    assert reviewed["assessment"]["step_assessments"][0]["assessment"] == assessed["assessment"]
    _, validity = call(workspace, "assess_retro_validity", **proposal["saved_candidate"])
    assert validity["proposal"]["mapped_reaction_smiles"] == candidate["condition_query_reaction_smiles"]
    assert validity["mapping_evidence"] == assessed["mapping_evidence"]
    _, revised = call(workspace, "revise_route_branch", source_ref=route.artifact_ref,
                       remove_step_ids=["s1"], replacement_steps=[{"external_step_id": "s2", **proposal}],
                       reason="Retain selected reconstruction", risks=["Experimental feasibility is unresolved"])
    assert revised["mapping_evidence"] == {"s2": assessed["mapping_evidence"]}
    assert workspace.store.read_artifact(workspace.replay(step.artifact_ref).artifact_ref)["matches"]


def _saved(workspace, *, mapping=None):
    candidate = {"realization_id": "real", "target_smiles": "CCN", "precursor_smiles": "CCBr.N",
                 "condition_query_reaction_smiles": mapping,
                 "forward_validation_status": "verified_signature"}
    event = workspace.store.append("call", {"operation": "disconnect_target", "execution_status": "completed",
        "result": {"strategies": [{"strategy_id": "strategy", "representative": candidate}]}})
    return {"target_smiles": "CCN", "precursor_smiles": "CCBr.N", "saved_candidate": {
        "source_ref": event.artifact_ref, "realization_id": "real", "strategy_id": "strategy"}}


@pytest.mark.parametrize("change,message", [
    ({"target_smiles": "CCO"}, "does not match target_smiles"),
    ({"precursor_smiles": "CCCl.N"}, "No saved candidate"),
    ({"mapped_reaction_smiles": "CCBr.N>>CCN"}, "not both"),
])
def test_candidate_conflicts_are_rejected_without_overwriting_evidence(workspace, change, message):
    proposal = _saved(workspace, mapping="[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]")
    event = workspace.run("assess_route_step", {"proposal": {**proposal, **change}})
    saved = workspace.store.read_artifact(event.artifact_ref)
    assert saved["execution_status"] == "error"
    assert message in saved["error"]["message"]


def test_candidate_status_never_bypasses_map_validation(workspace):
    proposal = _saved(workspace, mapping="[CH3:1][CH2:1][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]")
    _, result = call(workspace, "assess_route_step", proposal=proposal)
    assert result["assessment"]["admission_eligible"] is False
    assert any(gate["status"] not in {"pass", "not_run"} for gate in result["assessment"]["gates"])


def test_logged_epoxidation_reuses_saved_correspondence_without_claiming_feasibility(workspace):
    precursor = "C1=CCN2Cc3ccccc3[C@H]2C1.O=C(OO)c1cccc(Cl)c1"
    target = "c1ccc2c(c1)CN1CC3OC3C[C@H]21"
    mapping = ("[OH:11][O:902][C:901](=[O:900])[c:903]1[cH:904][cH:905][cH:906][c:907]([Cl:908])[cH:909]1."
               "[cH:1]1[cH:2][cH:3][c:4]2[c:5]([cH:6]1)[CH2:7][N:8]1[CH2:9][CH:10]=[CH:12][CH2:13][C@H:14]21>>"
               "[cH:1]1[cH:2][cH:3][c:4]2[c:5]([cH:6]1)[CH2:7][N:8]1[CH2:9][CH:10]3[O:11][CH:12]3[CH2:13][C@H:14]21")
    source = workspace.store.append("call", {"operation": "disconnect_target", "execution_status": "completed",
        "result": {"strategies": [{"strategy_id": "epoxidation", "representative": {
            "realization_id": "epoxide", "target_smiles": target, "precursor_smiles": precursor,
            "condition_query_reaction_smiles": mapping}}]}})
    proposal = {"target_smiles": target, "precursor_smiles": precursor}
    _, without = call(workspace, "assess_route_step", proposal=proposal)
    _, with_mapping = call(workspace, "assess_route_step", proposal={**proposal, "saved_candidate": {
        "source_ref": source.artifact_ref, "realization_id": "epoxide", "strategy_id": "epoxidation"}})
    def correspondence(result):
        return next(g for g in result["assessment"]["gates"] if g["gate_id"] == "atom_correspondence")["status"]
    assert correspondence(without) == "unresolved"
    assert correspondence(with_mapping) == "pass"
    assert with_mapping["experimental_feasibility"] == "not_established"


def test_shared_realization_needs_unambiguous_selection_and_mapping(workspace):
    payload = ambiguous_disconnection()
    event = workspace.store.append("call", payload)
    first = payload["result"]["strategies"][0]["representative"]
    from chem_coworker.scientific_workspace.adapters.step_selection import resolve_saved_candidate

    with pytest.raises(ValueError, match="Ambiguous"):
        resolve_saved_candidate(workspace.store, {"target_smiles": first["target_smiles"], "saved_candidate": {
            "source_ref": event.artifact_ref, "realization_id": "REAL2:shared"}})
    with pytest.raises(ValueError, match="no mapped reconstruction"):
        resolve_saved_candidate(workspace.store, {"target_smiles": first["target_smiles"],
            "precursor_smiles": first["precursor_smiles"], "saved_candidate": {
                "source_ref": event.artifact_ref, "realization_id": "REAL2:shared", "strategy_id": "s0"}})


def test_duplicate_candidate_with_conflicting_mapping_is_not_silently_selected(workspace):
    from chem_coworker.scientific_workspace.adapters.step_selection import resolve_saved_candidate

    candidate = {"realization_id": "real", "target_smiles": "CCN", "precursor_smiles": "CCBr.N",
                 "condition_query_reaction_smiles": "[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"}
    alternate = {**candidate, "condition_query_reaction_smiles": "[CH3:1][CH2:2][Br:3].[NH3:4]>>[CH3:2][CH2:1][NH2:4]"}
    event = workspace.store.append("call", {"operation": "disconnect_target", "execution_status": "completed",
        "result": {"strategies": [{"strategy_id": "strategy", "representative": candidate,
                                    "alternate_realizations": [alternate]}]}})
    with pytest.raises(ValueError, match="conflicting saved mappings"):
        resolve_saved_candidate(workspace.store, {"target_smiles": "CCN", "precursor_smiles": "CCBr.N",
            "saved_candidate": {"source_ref": event.artifact_ref, "realization_id": "real", "strategy_id": "strategy"}})
