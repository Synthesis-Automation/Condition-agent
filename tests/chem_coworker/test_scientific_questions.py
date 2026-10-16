"""Question provenance and saved dossier regressions, without live scientific providers."""

from copy import deepcopy

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.answers.answer_contracts import (
    ScientificAnswer, validate_answer_evidence,
)
from chem_coworker.scientific_workspace.console import ScientificConsole
from chem_coworker.scientific_workspace.core.store import canonical_bytes


@pytest.fixture
def workspace(tmp_path):
    store = InvestigationStore.create(tmp_path / "study", objective="Investigate a step", baseline={})
    return ScientificWorkspace(store.root)


def assessment(workspace, *, precursors="CCBr.N", status="ambiguous", step_id=None):
    """Save an assessment-shaped fixture; this does not test chemical correctness."""
    proposal = {"precursor_smiles": precursors, "target_smiles": "CCN"}
    checks = {"status": status, "admission_eligible": False, "actionable": False,
              "gates": [{"gate_id": "selectivity", "status": "not_run"}], "warnings": []}
    result = {"proposal": proposal, "assessment": checks}
    if step_id:
        result = {"proposal": {"steps": [{"external_step_id": step_id, **proposal}]},
                  "assessment": {"step_assessments": [{"external_step_id": step_id, "assessment": checks}]}}
    return workspace.store.append("call", {
        "operation": "assess_route_proposal" if step_id else "assess_route_step",
        "execution_status": "completed", "result": result,
    }).artifact_ref


def question(workspace, source, **kwargs):
    return workspace.record_question(
        source_ref=source, question="Is the competing site addressed?",
        resolution_criterion="Inspected evidence addresses the competing site",
        next_check="Inspect the nearest source experiment", **kwargs,
    )


def test_lifecycle_reopens_without_erasing_conflicts_or_old_resolution(workspace):
    source = assessment(workspace)
    q = question(workspace, source)
    conflict = assessment(workspace, status="conflicting")
    working = workspace.update_question(q.artifact_ref, status="investigating", rationale="Compare the conflict",
                                       contradicting_refs=(conflict,))
    resolved = workspace.update_question(working.artifact_ref, status="resolved",
                                        rationale="Criterion addressed within the stated scope; conflict remains a caution",
                                        resolution_refs=(source,))
    reopened = ScientificWorkspace(workspace.store.root)
    row = reopened.list_questions()["items"][0]
    assert row["status"] == "resolved" and row["review_status"] == "unreviewed"
    assert row["contradicting_refs"] == [conflict]
    old = canonical_bytes(reopened.store.read_artifact(resolved.artifact_ref))
    reopened.update_question(resolved.artifact_ref, status="open", rationale="New evidence requires another check")
    assert reopened.list_questions()["items"][0]["resolution_refs"] == []
    assert canonical_bytes(reopened.store.read_artifact(resolved.artifact_ref)) == old
    assert reopened.store.read_artifact(source)["result"]["assessment"]["admission_eligible"] is False


def test_duplicate_and_stale_updates_do_not_fork_history(workspace):
    q = question(workspace, assessment(workspace))
    with pytest.raises(ValueError, match="already exists"):
        question(workspace, workspace.store.read_artifact(q.artifact_ref)["source_ref"])
    workspace.update_question(q.artifact_ref, status="investigating", rationale="Read the source")
    count = len(workspace.store.events())
    with pytest.raises(ValueError, match="current question"):
        workspace.update_question(q.artifact_ref, status="withdrawn", rationale="A stale decision")
    assert len(workspace.store.events()) == count
    assert not (workspace.store.root / ".writer.lock").exists()


@pytest.mark.parametrize("kind,execution", [("decision", None), ("call", "timed_out"),
                                            ("call", "error"), ("call", "cancelled")])
def test_assertions_and_incomplete_calls_cannot_resolve_questions(workspace, kind, execution):
    q = question(workspace, assessment(workspace))
    ref = workspace.store.append(kind, {"execution_status": execution}).artifact_ref
    with pytest.raises(ValueError, match="completed scientific"):
        workspace.update_question(q.artifact_ref, status="resolved", rationale="Cannot prove the criterion",
                                  resolution_refs=(ref,))
    assert workspace.list_questions()["items"][0]["status"] == "open"


@pytest.mark.parametrize("changes", [
    {"status": "resolved", "rationale": "No evidence"},
    {"status": "withdrawn", "rationale": "   "},
    {"status": "withdrawn", "rationale": None},
    {"status": "verified", "rationale": "Not an allowed scientific claim"},
])
def test_invalid_updates_are_not_published(workspace, changes):
    q = question(workspace, assessment(workspace))
    count = len(workspace.store.events())
    with pytest.raises(ValueError):
        workspace.update_question(q.artifact_ref, **changes)
    assert len(workspace.store.events()) == count


def test_captured_passage_can_support_an_unreviewed_resolution(workspace):
    q = question(workspace, assessment(workspace))
    passage = workspace.capture_source("The isolated example retained the other functional group.",
                                       url="https://example.org/paper", locator="Example 1")
    workspace.update_question(q.artifact_ref, status="resolved", rationale="Source addresses the narrow criterion",
                              resolution_refs=(passage.artifact_ref,))
    assert workspace.list_questions()["resolution_authority"] == "agent_assessed_unreviewed_not_chemical_validation"


def test_exact_scope_and_pagination_do_not_transfer_resolutions(workspace):
    source = assessment(workspace, step_id="s1")
    q = question(workspace, source, step_id="s1")
    workspace.update_question(q.artifact_ref, status="resolved", rationale="An agent assessment",
                              resolution_refs=(source,))
    revised = assessment(workspace, precursors="CCCl.N", step_id="s1")
    assert workspace.list_questions(source_ref=revised, step_id="s1")["total"] == 0
    assert workspace.list_questions(source_ref=source)["total"] == 0
    with pytest.raises(ValueError, match="not present"):
        question(workspace, source, step_id="missing")
    question(workspace, revised, step_id="s1")
    page = workspace.list_questions(limit=1)
    assert page["total"] == 2 and page["next_offset"] == 1
    assert workspace.list_questions(offset=1, limit=1)["next_offset"] is None


def test_missing_resolution_evidence_fails_closed_on_reopen(workspace):
    source = assessment(workspace)
    q = question(workspace, source)
    evidence = assessment(workspace, status="conflicting")
    workspace.update_question(q.artifact_ref, status="resolved", rationale="Review the evidence", resolution_refs=(evidence,))
    (workspace.store.root / "artifacts" / (evidence.split(":")[1] + ".json")).unlink()
    with pytest.raises(OSError):
        ScientificWorkspace(workspace.store.root).list_questions()


def test_dossier_preserves_conflicts_exact_recipe_scope_and_history(workspace, monkeypatch):
    source = assessment(workspace)
    conflict = assessment(workspace, status="conflicting")
    q = question(workspace, source)
    def recipe(reaction):
        return workspace.store.append("call", {
            "operation": "assess_recipe", "execution_status": "completed",
            "arguments": {"reaction_smiles": reaction, "recipe": {"recipe_id": "recipe-A"}},
            "result": {"status": "conflicting", "hard_conflicts": ["fixture conflict"], "coverage": None},
        }).artifact_ref
    matching, unrelated = recipe("N.CCBr>>CCN"), recipe("CCCl.N>>CCN")
    monkeypatch.setattr(workspace.operations, "invoke", lambda *args: pytest.fail("Dossier must not run science"))
    event = workspace.step_dossier(source, evidence_refs=(matching, unrelated))
    dossier = workspace.store.read_artifact(event.artifact_ref)
    assert {a["status"] for a in dossier["assessments"]["structural_assessments"]} == {"ambiguous", "conflicting"}
    assert [a["artifact_ref"] for a in dossier["assessments"]["recipe_assessments"]] == [matching]
    assert dossier["experimental_feasibility"] == "not_established"
    assert conflict in event.evidence_refs and q.artifact_ref in event.evidence_refs
    original = deepcopy(dossier)
    workspace.update_question(q.artifact_ref, status="withdrawn", rationale="Superseded investigation")
    assert workspace.store.read_artifact(event.artifact_ref) == original
    assert workspace.call_summary(event)["inspection"]["result_path"] == "$"
    assert len(canonical_bytes(workspace.call_summary(event))) <= 16000


def test_recipe_dossier_retains_missing_checks(workspace):
    ref = workspace.store.append("call", {
        "operation": "assess_recipe", "execution_status": "completed",
        "arguments": {"reaction_smiles": "CCBr.N>>CCN", "recipe": {}},
        "result": {"status": "no_known_conflict", "compatible": True},
    }).artifact_ref
    dossier = workspace.store.read_artifact(workspace.step_dossier(ref).artifact_ref)
    assert dossier["assessments"]["structural_assessments"] == []
    assert dossier["assessments"]["recipe_assessments"][0]["coverage"] is None
    assert dossier["questions"]["total"] == 0


def test_console_discovery_and_no_question_or_dossier_as_scientific_citation(workspace):
    source = assessment(workspace)
    console = ScientificConsole(workspace)
    assert len(workspace.help(["record_question", "update_question", "list_questions", "step_dossier"])) == 4
    q = question(workspace, source)
    assert console.dispatch({"action": "list_questions"})["total"] == 1
    summary = console.dispatch({"action": "step_dossier", "arguments": {"source_ref": source}})
    for ref in (q.artifact_ref, summary["event"]["artifact_ref"]):
        draft = workspace.answer_template("Unreviewed investigation")
        draft["evidence_refs"] = [ref]
        with pytest.raises(ValueError, match="not agent assertions"):
            validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)
