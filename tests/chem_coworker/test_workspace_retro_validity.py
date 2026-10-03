"""Recorded validity discovery, source selection, bounded audits and replay."""

from dataclasses import asdict

import pytest

from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.shared_core_index import build_shared_core_index
from condition_recommender.sqlite_indexing import save_sqlite_generic_index
from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import artifact_identity
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from core_retrosynthesis import assess_retro_validity
from core_retrosynthesis.external_proposal_assessment import ExternalRetrosynthesisProposal
from tests.chem_coworker.test_route_investigation_workspace import library, workspace, call
from tests.condition_recommender.test_shared_core_retrieval import record
from tests.core_retrosynthesis_tests.test_external_proposal_admission import _route_value


PROPOSAL = {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"}


def test_tool_is_discoverable_and_records_canonical_axes(workspace, library):
    definition = workspace.operations.definition("assess_retro_validity")
    assert definition.required_artifacts == ("retro_library",)
    event, result = call(workspace, "assess_retro_validity", proposal=PROPOSAL)
    direct = assess_retro_validity(ExternalRetrosynthesisProposal.from_dict(PROPOSAL), library)
    assert canonical_bytes(result["validity"]) == canonical_bytes(direct.to_dict())
    summary = workspace.call_summary(event, detailed=True)
    assert summary["result_summary"]["validity"]["precedent_grade"] == "exact_reaction"
    brief = workspace.call_summary(event)
    assert brief["result_summary"]["validity"]["evidence_rank"] == 4
    replay = ScientificWorkspace(workspace.store.root).replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"]


def test_saved_single_step_and_validity_sources_are_linked_and_inspectable(workspace):
    source, _ = call(workspace, "assess_route_step", proposal=PROPOSAL)
    assessed, result = call(workspace, "assess_retro_validity", source_ref=source.artifact_ref)
    assert assessed.evidence_refs == (source.artifact_ref,)
    assert result["validity"]["evidence_rank"] == 4
    _, inspected = call(workspace, "inspect_step_precedents", source_ref=assessed.artifact_ref)
    assert inspected["precedents"]
    _, again = call(workspace, "assess_retro_validity", source_ref=assessed.artifact_ref)
    assert again["validity"] == result["validity"]


def test_saved_disconnection_and_route_selection_bind_concrete_realizations(workspace):
    source, disconnected = call(workspace, "disconnect_target", target_smiles="CCN", top_k=2,
                                max_templates_to_apply=10, max_candidates_to_validate=5)
    candidate = disconnected["strategies"][0]["representative"]
    _, result = call(workspace, "assess_retro_validity", source_ref=source.artifact_ref,
                     realization_id=candidate["realization_id"])
    assert result["validity"]["structural_status"] == "verified"
    route, _ = call(workspace, "assess_route_proposal", proposal=_route_value())
    _, step = call(workspace, "assess_retro_validity", source_ref=route.artifact_ref, step_id="step-1")
    assert step["validity"]["precedent_grade"] == "exact_reaction"


def test_optional_pinned_corpus_is_used_without_recipe_ranking(workspace, tmp_path, monkeypatch):
    index = build_generic_index([record(1, "CC=O.N>>CCN")])
    index_path, shared_path = tmp_path / "index.sqlite", tmp_path / "shared.sqlite"
    save_sqlite_generic_index(index, index_path)
    from condition_recommender.sqlite_indexing import load_sqlite_generic_index

    index = load_sqlite_generic_index(index_path)
    build_shared_core_index(index, shared_path)
    # Construct a fresh investigation rather than modifying the running baseline.
    from chem_coworker.scientific_workspace import InvestigationStore

    baseline = {**workspace.store.manifest["baseline"], "artifacts": {
        **workspace.store.manifest["baseline"]["artifacts"],
        "condition_index": artifact_identity(index_path), "shared_core_index": artifact_identity(shared_path),
    }}
    store = InvestigationStore.create(tmp_path / "corpus-run", objective="Grade a precedent", baseline=baseline)
    corpus_workspace = ScientificWorkspace(store.root)
    from condition_recommender import GenericConditionRecommender

    monkeypatch.setattr(GenericConditionRecommender, "recommend", lambda *args, **kwargs:
                        pytest.fail("Evidence grading must not run recipe ranking"))
    _, result = call(corpus_workspace, "assess_retro_validity", proposal=PROPOSAL)
    support = result["validity"]["corpus_precedent_support"]
    assert support["strongest_level"] == "whole_reaction"
    assert support["matches"][0]["yield_pct"] == 70


def saved_forward(workspace, route_ref, *, execution="completed", target="CCN"):
    from core_retrosynthesis import RetroForwardEvidence

    projection = asdict(RetroForwardEvidence(
        starting_materials="CC=O.N", intended_product=target,
        validity="structurally_supported_with_competition", targeted_replay_status="structurally_reproduced",
        intended_match="exact", best_competitor_product="CCO", warnings=("COMPETING_FORWARD_PRODUCTS",),
    ))
    return workspace.store.append("call", {
        "operation": "assess_route_step_forward", "execution_status": execution,
        "result": {"source_ref": route_ref, "step_id": "step-1", "execution_status": execution,
                   "assessment": projection if execution == "completed" else None,
                   "error": None if execution == "completed" else {"type": "TimeoutExpired"}},
    })


@pytest.mark.parametrize("execution", ["completed", "timed_out", "error", "cancelled"])
def test_saved_forward_audit_preserves_competition_or_partial_execution(workspace, execution):
    route, _ = call(workspace, "assess_route_proposal", proposal=_route_value())
    audit = saved_forward(workspace, route.artifact_ref, execution=execution)
    event, result = call(workspace, "assess_retro_validity", proposal=PROPOSAL, forward_ref=audit.artifact_ref)
    assert audit.artifact_ref in event.evidence_refs
    assert result["validity"]["forward_execution_status"] == execution
    assert result["validity"]["evidence_rank"] == 4
    if execution == "completed":
        assert result["validity"]["status"] == "supported_with_cautions"
    else:
        assert f"forward:{execution}" in result["validity"]["unresolved_checks"]
        assert result["forward_execution"]["error"]["type"] == "TimeoutExpired"


def test_forward_evidence_for_other_structure_is_a_recorded_error(workspace):
    route, _ = call(workspace, "assess_route_proposal", proposal=_route_value())
    audit = saved_forward(workspace, route.artifact_ref, target="CCO")
    event = workspace.run("assess_retro_validity", {"proposal": PROPOSAL, "forward_ref": audit.artifact_ref})
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "error"
    assert "match this realization" in record["error"]["message"]


@pytest.mark.parametrize("arguments", [
    {}, {"proposal": PROPOSAL, "source_ref": "sha256:missing"},
    {"proposal": PROPOSAL, "realization_id": "wrong"},
    {"proposal": PROPOSAL, "candidate_limit": True},
])
def test_invalid_selections_and_limits_are_recorded_failures(workspace, arguments):
    event = workspace.run("assess_retro_validity", arguments)
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "error"
