"""Source attribution, exact-draft self-review, and environment diagnostics."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.capabilities import local_capabilities
from chem_coworker.scientific_workspace.evidence_review import REVIEW_AREAS, matching_evidence_review
from chem_coworker.scientific_workspace.store import InvestigationStore


@pytest.fixture
def research(tmp_path: Path):
    store = InvestigationStore.create(tmp_path / "research", objective="Evidence regression", baseline={})
    workspace = ScientificWorkspace(store.root)
    source = workspace.capture_source(
        "Example 2: the cis racemate was isolated in 80% yield. No enantiomeric purity reported.",
        url="https://example.org/patent", title="Synthetic test fixture", locator="Example 2",
    )
    excerpt = workspace.record_source_excerpt(source.artifact_ref, excerpt="the cis racemate", locator="Example 2")
    answer = ScientificAnswer.model_validate({
        "schema_version": "scientific_answer.v2",
        "answer_markdown": "A racemate precedent does not establish a synthesis of one enantiomer.",
        "evidence_refs": [], "uncertainties": ["Enantiomeric purity unknown"], "needs_user_input": False,
        "sources": [{"id": "p1", "kind": "external_source", "title": "Synthetic test fixture",
                     "artifact_ref": excerpt.artifact_ref, "url": "https://example.org/patent#example2",
                     "locator": "Example 2"}],
        "molecules": [], "target_molecule_ids": [], "steps": [], "routes": [],
        "claims": [{"text": "The supplied passage describes a racemate", "basis": "reported",
                    "source_ids": ["p1"], "limitations": ["Agent-supplied test excerpt"]}],
    })
    return workspace, source, excerpt, answer


def findings(reference: str) -> list[dict]:
    return [{
        "area": area, "claim": "The exact enantiomer is not established by this racemate passage",
        "assessment": "partial" if area == "structure_and_stereochemistry" else "not_checked",
        "evidence_refs": [reference] if area == "structure_and_stereochemistry" else [],
        "reason": "Synthetic regression: no independent source or experiment checked",
    } for area in REVIEW_AREAS]


def test_workspace_source_tools_and_precise_url_attribution(research) -> None:
    workspace, source, excerpt, answer = research
    view = workspace.inspect_source(source.artifact_ref, query="racemate", limit=100)
    assert view["query_found"] and "cis racemate" in view["text"]
    assert view["acquisition"] == "agent_supplied_excerpt"
    assert validate_answer_evidence(answer, workspace.store) == [excerpt.artifact_ref]
    answer.sources[0].url = "https://example.org/a-different-patent"
    with pytest.raises(ValueError, match="does not match"):
        validate_answer_evidence(answer, workspace.store)


def test_failed_retrieval_is_citable_failure_but_not_a_reported_passage(research) -> None:
    workspace, _, _, answer = research
    failure = workspace.store.append("literature_source", {
        "source_url": "https://example.org/patent", "retrieval_status": "failed",
        "extraction": {"status": "not_attempted", "text": ""},
    })
    answer.sources[0].artifact_ref = failure.artifact_ref
    with pytest.raises(ValueError, match="captured text"):
        validate_answer_evidence(answer, workspace.store)
    answer.sources = []
    answer.claims = []
    answer.evidence_refs = [failure.artifact_ref]
    assert validate_answer_evidence(answer, workspace.store) == [failure.artifact_ref]


def test_literature_cannot_be_relabelled_as_a_computation(research) -> None:
    workspace, _, _, answer = research
    answer.sources[0].kind = "local_artifact"
    answer.sources[0].url = None
    answer.claims[0].basis = "computed"
    with pytest.raises(ValueError, match="completed recorded"):
        validate_answer_evidence(answer, workspace.store)


def test_self_review_is_bound_to_exact_answer_and_does_not_certify_chemistry(research) -> None:
    workspace, _, excerpt, answer = research
    review = workspace.record_evidence_review(answer.model_dump(), findings(excerpt.artifact_ref))
    saved = workspace.store.read_artifact(review.artifact_ref)
    assert saved["origin"] == "agent_authored_self_review"
    assert saved["semantic_verification"] == "not_performed_by_validator"
    assert saved["unresolved_count"] == 5
    answer.evidence_refs = validate_answer_evidence(answer, workspace.store)
    assert matching_evidence_review(workspace.store, answer) == review.artifact_ref
    next_turn = workspace.store.append("user_message", {"text": "A new constraint needs another review"})
    assert matching_evidence_review(workspace.store, answer, after_sequence=next_turn.sequence) is None
    answer.answer_markdown += " A subsequently changed recommendation."
    assert matching_evidence_review(workspace.store, answer) is None
    answer.evidence_refs = [review.artifact_ref]
    with pytest.raises(ValueError, match="scientific evidence"):
        validate_answer_evidence(answer, workspace.store)


def test_self_review_requires_explicit_missing_checks_and_real_evidence(research) -> None:
    workspace, _, excerpt, answer = research
    checks = findings(excerpt.artifact_ref)
    with pytest.raises(ValueError, match="every area"):
        workspace.record_evidence_review(answer.model_dump(), checks[:-1])
    unsupported = deepcopy(checks)
    unsupported[0]["assessment"] = "supported"
    with pytest.raises(ValueError, match="require evidence_refs"):
        workspace.record_evidence_review(answer.model_dump(), unsupported)
    note = workspace.store.note("review", "Author assertion without source")
    checks[1]["evidence_refs"] = [note.artifact_ref]
    with pytest.raises(ValueError, match="not agent assertions"):
        workspace.record_evidence_review(answer.model_dump(), checks)


def test_capability_probe_does_not_promote_file_presence_or_requested_web_access(tmp_path: Path) -> None:
    present = tmp_path / "not-an-index.sqlite"
    present.write_text("This is deliberately not a database", "utf-8")
    report = local_capabilities({"artifacts": {
        "index": {"path": str(present), "status": "present"},
        "missing": {"path": str(tmp_path / "missing"), "status": "missing"},
    }})
    assert report["rdkit"]["status"] == "available"
    assert report["artifacts"]["index"]["status"] == "file_present"
    assert report["artifacts"]["index"]["content_validation"] == "not_checked"
    assert report["artifacts"]["missing"]["status"] == "missing"
    assert report["agent_web_search"]["status"] == "not_checked"
    assert report["model_and_reasoning"]["status"] == "not_confirmed"
