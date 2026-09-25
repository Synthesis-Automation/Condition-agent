"""Recorded proposed-route assessment, branch revision, recovery and comparison."""

from copy import deepcopy
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.baseline import artifact_identity, code_manifest, environment_versions
from chem_coworker.scientific_workspace.store import canonical_bytes
from core_retrosynthesis.external_proposal_assessment import ExternalRetrosynthesisProposal, assess_external_retrosynthesis_proposal
from core_retrosynthesis.external_route_admission import ExternalRouteProposal, assess_external_route_proposal
from core_retrosynthesis.generic_library import build_generic_library, save_generic_library
from tests.core_retrosynthesis_tests.test_external_proposal_admission import _row, _route_value, FIRST_REACTION, SECOND_REACTION


@pytest.fixture(scope="module")
def library():
    return build_generic_library(tuple(_row(reaction, ordinal) for ordinal, reaction in enumerate(
        (FIRST_REACTION, SECOND_REACTION, "CCBr.N>>CCN"), 1,
    )), levels=("L0", "L1", "L2"))


@pytest.fixture
def workspace(tmp_path: Path, library) -> ScientificWorkspace:
    root = Path(__file__).resolve().parents[2]
    path = tmp_path / "operators.json.gz"
    save_generic_library(library, path)
    baseline = {"repository": str(root), "code_files": code_manifest(root), "environment": environment_versions(),
                "artifacts": {"retro_library": artifact_identity(path)}}
    InvestigationStore.create(tmp_path / "investigation", objective="Development route revision", baseline=baseline)
    return ScientificWorkspace(tmp_path / "investigation")


def call(workspace, operation, **arguments):
    event = workspace.run(operation, arguments)
    saved = workspace.store.read_artifact(event.artifact_ref)
    assert saved["execution_status"] == "completed", saved
    return event, saved["result"]


@pytest.mark.parametrize("proposal", [
    {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"},
    {"target_smiles": "CCN", "precursor_smiles": "N.CC=O"},
    {"target_smiles": "CCN", "precursor_smiles": "CC.CN"},
    {"target_smiles": "invalid", "precursor_smiles": "CCO"},
    {"target_smiles": "CCO", "precursor_smiles": "CC=O.N",
     "mapped_reaction_smiles": "[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"},
])
def test_step_adapter_preserves_positive_negative_ambiguous_and_conflicting_evidence(workspace, library, proposal) -> None:
    _, result = call(workspace, "assess_route_step", proposal=proposal)
    direct = assess_external_retrosynthesis_proposal(ExternalRetrosynthesisProposal.from_dict(proposal), library)
    assert canonical_bytes(result["assessment"]) == canonical_bytes(direct.to_dict())
    assert result["experimental_feasibility"] == "not_established"
    assert result["assessment"]["actionable"] is False


def test_route_assessment_and_inspection_preserve_domain_evidence(workspace, library) -> None:
    event, record = call(workspace, "assess_route_proposal", proposal=_route_value())
    direct = assess_external_route_proposal(ExternalRouteProposal.from_dict(_route_value()), library)
    assert canonical_bytes(record["assessment"]) == canonical_bytes(direct.to_dict())
    inspected, view = call(workspace, "inspect_route_step", source_ref=event.artifact_ref, step_id="step-1")
    assert view["downstream_step_ids"] == ["step-2"]
    assert "CCN" in view["molecule_audits"]
    assert view["assessment"]["precedent_matches"]
    assert inspected.evidence_refs == (event.artifact_ref,)


def test_branch_revision_rechecks_preserved_downstream_step_and_replays_after_reopen(workspace, monkeypatch) -> None:
    source, before = call(workspace, "assess_route_proposal", proposal=_route_value(), unavailable_starting_materials=["CC=O"])
    source_bytes = canonical_bytes(workspace.store.read_artifact(source.artifact_ref))
    from core_retrosynthesis import external_route_admission

    original = external_route_admission.assess_external_retrosynthesis_proposal
    checked = []

    def spy(proposal, *args, **kwargs):
        checked.append(proposal.target_smiles)
        return original(proposal, *args, **kwargs)

    monkeypatch.setattr(external_route_admission, "assess_external_retrosynthesis_proposal", spy)
    reopened = ScientificWorkspace(workspace.store.root)
    revision, after = call(reopened, "revise_route_branch", source_ref=source.artifact_ref,
                          remove_step_ids=["step-1"], replacement_steps=[{
                              "external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CCBr.N",
                          }], reason="Declared aldehyde is unavailable; investigate an alternative amine preparation",
                          risks=["Overalkylation and experimental feasibility remain unverified"])
    assert len(checked) == 2 and _route_value()["target_smiles"] in checked
    assert before["material_constraints"]["status"] == "violated"
    assert after["material_constraints"]["status"] == "satisfied_for_declared_constraints"
    assert after["revision"]["preserved_step_ids"] == ["step-2"]
    assert after["revision"]["improvement_status"] == "not_automatically_established"
    assert after["assessment_options"] == before["assessment_options"]
    assert canonical_bytes(workspace.store.read_artifact(source.artifact_ref)) == source_bytes
    assert source.artifact_ref in revision.evidence_refs
    compare, table = call(reopened, "compare_route_proposals", source_refs=[source.artifact_ref, revision.artifact_ref])
    assert len(table["alternatives"]) == 2 and table["ranking"] == "not_performed"
    assert compare.evidence_refs == (source.artifact_ref, revision.artifact_ref)
    replay = ScientificWorkspace(workspace.store.root).replay(revision.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_unsupported_revision_remains_inspectable_and_noop_is_an_error(workspace) -> None:
    source, _ = call(workspace, "assess_route_proposal", proposal=_route_value())
    unknown, record = call(workspace, "revise_route_branch", source_ref=source.artifact_ref,
                          remove_step_ids=["step-1"], replacement_steps=[{
                              "external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CC.CN",
                          }], reason="Investigate an unsupported alternative", risks=["No validated correspondence"])
    assert record["assessment"]["status"] == "partially_supported"
    assert "step-1" in record["assessment"]["unresolved_step_ids"]
    _, broken = call(workspace, "revise_route_branch", source_ref=unknown.artifact_ref,
                     remove_step_ids=["step-2"], replacement_steps=[], reason="Demonstrate disconnected topology",
                     risks=["Target no longer produced"])
    assert broken["assessment"]["status"] == "invalid"
    assert broken["material_constraints"]["status"] == "unknown"
    error = workspace.run("revise_route_branch", {
        "source_ref": source.artifact_ref, "remove_step_ids": [], "replacement_steps": [],
        "reason": "No change", "risks": ["No change"],
    })
    assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"


def test_unknown_fields_bad_sources_and_mismatched_comparisons_are_rejected(workspace) -> None:
    value = _route_value()
    altered = deepcopy(value)
    altered["steps"][0]["confidence"] = 1
    for proposal in (altered, {**value, "admission_status": "verified"}):
        error = workspace.run("assess_route_proposal", {"proposal": proposal})
        assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"
    good, _ = call(workspace, "assess_route_proposal", proposal=value)
    constrained, _ = call(workspace, "assess_route_proposal", proposal=value, unavailable_starting_materials=["CCN"])
    failed = workspace.run("compare_route_proposals", {"source_refs": [good.artifact_ref, constrained.artifact_ref]})
    assert "same options" in workspace.store.read_artifact(failed.artifact_ref)["error"]["message"]
    note = workspace.store.note("hypothesis", "Not scientific execution evidence")
    failed = workspace.run("assess_route_proposal", {"proposal": value, "evidence_refs": [note.artifact_ref]})
    assert workspace.store.read_artifact(failed.artifact_ref)["execution_status"] == "error"


def test_supplied_recipe_is_assessed_separately_and_forward_challenge_is_optional(workspace) -> None:
    _, recipe = call(workspace, "resolve_recipe", components=[{
        "raw_identifier": "ethanol", "source_field": "user", "identifier_type": "name",
    }])
    _, result = call(workspace, "assess_route_step", proposal={
        "target_smiles": "CCN", "precursor_smiles": "CC=O.N", "proposed_conditions": recipe,
    }, include_forward=True)
    assert result["proposed_recipe_assessment"]["schema_version"] == "reaction_recipe_assessment.v2"
    assert result["assessment"]["forward_assessment"] is not None
    assert next(gate for gate in result["assessment"]["gates"] if gate["gate_id"] == "condition_support")["status"] == "unresolved"


def test_prepare_selected_serialized_tree_preserves_source_and_reassesses(workspace, library) -> None:
    from core_retrosynthesis.external_route_admission import external_route_proposal_from_tree

    # The planner and admission system share the same typed route-tree contract.
    tree = assess_external_route_proposal(ExternalRouteProposal.from_dict(_route_value()), library).admitted_route_tree
    assert tree is not None
    source = workspace.store.append("call", {
        "operation": "plan_routes", "execution_status": "completed", "arguments": {},
        "result": {"response": {"result": {"routes": [{"route_id": "selected", "route_tree": tree.to_dict()}], "partial_routes": []}}},
    })
    event, record = call(workspace, "prepare_route_proposal", source_ref=source.artifact_ref,
                         route_id="selected", unavailable_starting_materials=["CC=O"])
    direct = assess_external_route_proposal(external_route_proposal_from_tree(tree), library)
    assert canonical_bytes(record["assessment"]) == canonical_bytes(direct.to_dict())
    assert record["planner_route_id"] == "selected"
    assert record["material_constraints"]["status"] == "violated"
    assert event.evidence_refs == (source.artifact_ref,)
    failed = workspace.run("prepare_route_proposal", {"source_ref": source.artifact_ref, "route_id": "invented"})
    assert workspace.store.read_artifact(failed.artifact_ref)["execution_status"] == "error"


def test_material_constraint_cannot_be_silently_removed_during_revision(workspace) -> None:
    source, _ = call(workspace, "assess_route_proposal", proposal=_route_value(), unavailable_starting_materials=["CC=O"])
    error = workspace.run("revise_route_branch", {
        "source_ref": source.artifact_ref, "remove_step_ids": [], "replacement_steps": [],
        "reason": "Attempt to change constraints", "risks": ["Unresolved"], "unavailable_starting_materials": [],
    })
    assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"


def test_requested_condition_retrieval_requires_recorded_data(workspace) -> None:
    error = workspace.run("assess_route_step", {
        "proposal": {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"}, "include_conditions": True,
    })
    saved = workspace.store.read_artifact(error.artifact_ref)
    assert saved["execution_status"] == "error"
    assert "condition_index" in saved["error"]["message"]
