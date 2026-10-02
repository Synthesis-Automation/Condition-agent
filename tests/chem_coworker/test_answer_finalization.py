"""Compact authoring preserves attribution, validated handoff and exact self-review."""

from copy import deepcopy
import json

import pytest

from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.answers.answer_contracts import (
    ScientificAnswer,
    validate_answer_evidence,
)
from chem_coworker.scientific_workspace.answers.answer_handoff import load_answer_handoff
from chem_coworker.scientific_workspace.answers.evidence_review import (
    REVIEW_AREAS,
    matching_evidence_review,
)
from chem_coworker.scientific_workspace.core.store import InvestigationStore


@pytest.fixture
def workspace(tmp_path):
    root = tmp_path / "investigation"
    InvestigationStore.create(root, objective="Answer authoring test", baseline={})
    attempt = root / "turns" / "attempt-1"
    attempt.mkdir(parents=True)
    return ScientificWorkspace(root), attempt / "answer-draft.json"


def test_compact_draft_roundtrips_through_existing_handoff_and_exact_review(workspace):
    w, path = workspace
    source = w.capture_source("Example: conversion is reported, without a yield.",
                              url="https://example.org/paper", locator="Example 1")
    draft = {
        "answer_markdown": "A proposed conversion; its yield remains unknown.",
        "sources": [{"id": "paper", "kind": "external_source", "title": "Example 1",
                     "url": "https://example.org/paper", "locator": "Example 1",
                     "artifact_ref": source.artifact_ref}],
        "molecules": [{"id": "a", "name": "Input", "smiles": "CCO", "basis": "input"},
                      {"id": "b", "name": "Proposed product", "smiles": "CC=O", "basis": "proposed"}],
        "target_molecule_ids": ["b"],
        "steps": [{"id": "s1", "title": "Proposed conversion", "basis": "proposed",
                   "reactant_ids": ["a"], "product_ids": ["b"], "source_ids": ["paper"],
                   "conditions": [{"text": "Conditions unknown", "basis": "unknown"}]}],
        "routes": [{"id": "r1", "title": "Proposal", "step_ids": ["s1"]}],
    }
    before = deepcopy(draft)
    findings = [{"area": area, "claim": "Fixture is not chemically verified.",
                 "assessment": "not_checked", "evidence_refs": [],
                 "reason": "Only answer transport is under test."} for area in REVIEW_AREAS]
    receipt = w.finalize_answer(path, draft, findings=findings)
    (path.parent / "agent-final.json").write_text(json.dumps(receipt), encoding="utf-8")
    saved = load_answer_handoff(w.store.root, path.parent, thread_id="fixture", usage={})
    answer = ScientificAnswer.model_validate(saved)
    assert draft == before
    assert answer.steps[0].yield_info is None
    assert answer.steps[0].basis == "proposed"
    assert answer.evidence_refs == validate_answer_evidence(answer, w.store) == [source.artifact_ref]
    assert answer.claims == [] and answer.uncertainties == []
    review_ref = matching_evidence_review(w.store, answer)
    assert review_ref
    assert w.store.read_artifact(review_ref)["semantic_verification"] == "not_performed_by_validator"
    answer.answer_markdown = "Changed conclusion"
    assert matching_evidence_review(w.store, answer) is None


@pytest.mark.parametrize("change", [
    {"claims": [{"text": "Unsupported claim", "basis": "reported"}]},
    {"claims": [{"text": "Missing attribution"}]},
    {"uncertainties": None},
    {"schema_version": "wrong"},
    {"extra": "must not be discarded"},
    {"sources": [{"id": "paper", "kind": "external_source", "title": "Unlinked",
                  "artifact_ref": "sha256:" + "0" * 64, "locator": "Example 1"}]},
])
def test_defaults_do_not_relax_validation_or_overwrite_valid_draft(workspace, change):
    w, path = workspace
    w.finalize_answer(path, {"answer_markdown": "A clarification is needed.", "needs_user_input": True})
    original = path.read_bytes()
    with pytest.raises(ValueError):
        w.finalize_answer(path, {"answer_markdown": "Invalid draft", **change})
    assert path.read_bytes() == original
    assert w.store.events() == ()  # No fabricated computation or self-review.


def test_failed_source_and_bad_review_are_rejected_before_saving(workspace):
    w, path = workspace
    failed = w.store.append("call", {"operation": "analyze_molecule", "execution_status": "error"})
    draft = {"answer_markdown": "Computed claim", "claims": [
        {"text": "Claim", "basis": "computed", "source_ids": ["check"]}],
        "sources": [{"id": "check", "kind": "local_artifact", "title": "Failed check",
                     "artifact_ref": failed.artifact_ref, "locator": "result"}]}
    with pytest.raises(ValueError, match="completed recorded"):
        w.finalize_answer(path, draft)
    with pytest.raises(ValueError, match="explicit review findings"):
        w.finalize_answer(path, {"answer_markdown": "Unreviewed"}, findings=[])
    assert not path.exists()
    assert len(w.store.events()) == 1


@pytest.mark.parametrize("destination", ["outside", "evidence"])
def test_submission_cannot_write_outside_workspace_or_overwrite_artifacts(workspace, destination):
    w, path = workspace
    target = w.store.root.parent / "answer-draft.json" if destination == "outside" else w.store.root / "investigation.json"
    before = target.read_bytes() if target.exists() else None
    with pytest.raises(ValueError):
        w.finalize_answer(target, {"answer_markdown": "Must not be written"})
    assert (target.read_bytes() if target.exists() else None) == before
