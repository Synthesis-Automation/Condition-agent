"""Traceable agent self-review; structural checks never imply semantic verification."""

from __future__ import annotations

import hashlib
from typing import Any, Literal, Mapping

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .answer_contracts import ScientificAnswer, validate_answer_evidence
from .store import InvestigationEvent, InvestigationStore, canonical_bytes


REVIEW_AREAS = (
    "source_identity", "structure_and_stereochemistry", "conditions_and_yields",
    "route_completeness", "counterevidence",
)


class ReviewFinding(BaseModel):
    """An agent-authored assessment, not an automatically validated chemical finding."""

    model_config = ConfigDict(extra="forbid", strict=True)
    area: Literal[
        "source_identity", "structure_and_stereochemistry", "conditions_and_yields",
        "route_completeness", "counterevidence",
    ]
    claim: str = Field(min_length=1, max_length=4000)
    assessment: Literal["supported", "partial", "unsupported", "conflicting", "not_checked", "not_applicable"]
    evidence_refs: list[str] = Field(max_length=30)
    reason: str = Field(min_length=1, max_length=4000)

    @model_validator(mode="after")
    def require_evidence(self) -> "ReviewFinding":
        """Supporting or conflicting evidence assessments must name their artifacts."""
        if self.assessment in {"supported", "partial", "conflicting"} and not self.evidence_refs:
            raise ValueError("Evidence assessments require evidence_refs")
        return self


def answer_digest(draft: Mapping[str, Any]) -> str:
    """Identify the exact typed draft without depending on citation list ordering."""
    normalized = ScientificAnswer.model_validate(dict(draft)).model_dump()
    normalized["evidence_refs"] = sorted(set(normalized["evidence_refs"]))
    return hashlib.sha256(canonical_bytes(normalized)).hexdigest()


def record_evidence_review(
    store: InvestigationStore, draft: Mapping[str, Any], findings: list[dict[str, Any]],
) -> InvestigationEvent:
    """Save a complete checklist with missing work explicit; never promote its findings.

    The content hash binds the review to the draft. A subsequent answer revision
    requires a new review to count as reviewed by its author. No semantic claim
    is checked here; exact source membership belongs to literature excerpts.
    """
    if not 1 <= len(findings) <= 50:
        raise ValueError("Supply 1 to 50 explicit review findings")
    answer = ScientificAnswer.model_validate(dict(draft))
    answer.evidence_refs = validate_answer_evidence(answer, store)
    parsed = [ReviewFinding.model_validate(item) for item in findings]
    if {item.area for item in parsed} != set(REVIEW_AREAS):
        raise ValueError("Review must cover every area; use not_checked or not_applicable explicitly")
    allowed = {event.artifact_ref for event in store.events() if event.kind in {
        "call", "derived_file", "replay", "custom_execution", "literature_source", "literature_excerpt",
    }}
    evidence = set(answer.evidence_refs)
    for item in parsed:
        for reference in item.evidence_refs:
            store.read_artifact(reference)
            if reference not in allowed:
                raise ValueError("Review findings must link scientific evidence, not agent assertions")
            evidence.add(reference)
    return store.append("evidence_review", {
        "schema_version": "scientific_evidence_review.v1",
        "answer_sha256": answer_digest(answer.model_dump()),
        "draft": answer.model_dump(),
        "origin": "agent_authored_self_review", "review_status": "not_independently_reviewed",
        "semantic_verification": "not_performed_by_validator",
        "findings": [item.model_dump() for item in parsed],
        "unresolved_count": sum(item.assessment in {
            "partial", "unsupported", "conflicting", "not_checked",
        } for item in parsed),
    }, evidence_refs=tuple(sorted(evidence)))


def matching_evidence_review(
    store: InvestigationStore, answer: ScientificAnswer, *, after_sequence: int = 0,
) -> str | None:
    """Find a self-review for this exact final draft; stale reviews do not count."""
    digest = answer_digest(answer.model_dump())
    for event in reversed(store.events()):
        if event.sequence <= after_sequence:
            break
        if event.kind == "evidence_review":
            value = store.read_artifact(event.artifact_ref)
            if value.get("answer_sha256") == digest:
                return event.artifact_ref
    return None
