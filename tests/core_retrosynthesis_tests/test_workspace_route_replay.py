"""Persisted workspace route revisions retain authoritative domain checks."""

from pathlib import Path

from chem_coworker.contracts import MultistepRetrosynthesisRequest, MultistepRetrosynthesisResponse
from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.baseline import code_manifest, environment_versions
from chem_coworker.scientific_workspace.operations import ScientificOperations
from core_retrosynthesis import plan_multistep_routes

from .test_multistep import _LiteratureIndex, _expander
from .test_route_refinement import _condition_evaluator, _identified_candidate


def test_route_revision_rehydrates_source_after_session_change(tmp_path, monkeypatch) -> None:
    source_candidate = _identified_candidate(
        "CCCCCCCC", "C.C", "source", realization="realization:source", score=0.99,
    )
    alternative = _identified_candidate(
        "CCCCCCCC", "N.O", "alternative", realization="realization:alternative", score=0.90,
    )
    expander = _expander({"CCCCCCCC": (source_candidate, alternative)})
    invocations = []

    def plan(self, settings, exclusions=()):
        invocations.append(exclusions)
        result = plan_multistep_routes(
            "CCCCCCCC", object(), _LiteratureIndex(), max_depth=2,
            molecular_weight_threshold=50.0, top_k_routes=2,
            condition_evidence_evaluator=_condition_evaluator,
            candidate_exclusions=exclusions, expander=expander,
        )
        return MultistepRetrosynthesisResponse(
            request=MultistepRetrosynthesisRequest(**settings), valid=True, result=result,
        )

    monkeypatch.setattr(ScientificOperations, "_plan", plan)
    root = Path(__file__).resolve().parents[2]
    store = InvestigationStore.create(tmp_path / "run", objective="Test recorded refinement",
                                      baseline={"repository": str(root), "code_files": code_manifest(root),
                                                "environment": environment_versions(), "artifacts": {}})
    workspace = ScientificWorkspace(store.root)
    first = workspace.run("plan_routes", {"settings": {"target_smiles": "CCCCCCCC"}})
    original = store.read_artifact(first.artifact_ref)
    assert original["execution_status"] == "completed"
    assessment = next(item for item in original["result"]["assessments"]
                      if any(issue["kind"] == "condition_gap" for issue in item["issues"]))
    issue = next(item for item in assessment["issues"] if item["kind"] == "condition_gap")
    intent = {"source_route_id": assessment["route_id"], "source_step_id": issue["subject_id"],
              "objective": "resolve_condition_gap", "method": "alternate_disconnection",
              "issue_ids": [issue["issue_id"]]}
    reopened = ScientificWorkspace(store.root)
    revised = reopened.run("revise_routes", {"source_ref": first.artifact_ref, "intent": intent})
    result = store.read_artifact(revised.artifact_ref)
    assert result["execution_status"] == "completed"
    assert result["result"]["refinement"]["status"] == "improved_alternative_found"
    assert result["result"]["refinement"]["source_route_preserved"] is True
    assert store.read_artifact(first.artifact_ref) == original
    assert revised.evidence_refs == (first.artifact_ref,)
    assert len(invocations) == 3  # Initial run, rehydration, revised search.
    assert invocations[-1][0].strategy_id == source_candidate.strategy_id

    invalid = reopened.run("revise_routes", {
        "source_ref": first.artifact_ref, "intent": {**intent, "issue_ids": ["invented_issue"]},
    })
    assert store.read_artifact(invalid.artifact_ref)["execution_status"] == "error"
    assert "unknown issue" in store.read_artifact(invalid.artifact_ref)["error"]["message"]

    drifted = ScientificWorkspace(store.root)
    monkeypatch.setattr(drifted.operations, "invoke", lambda operation, arguments: {"changed": True})
    # Call the method directly so its internal rehydration uses the changed dispatcher.
    import pytest
    with pytest.raises(ValueError, match="replay differs"):
        drifted.operations.revise_routes(first.artifact_ref, intent)
