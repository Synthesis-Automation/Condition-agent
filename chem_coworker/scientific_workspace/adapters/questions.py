"""Evidence-linked, agent-assessed questions over immutable investigation records."""

from __future__ import annotations

import hashlib
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..core.store import InvestigationEvent, InvestigationStore, canonical_bytes
from .literature import load_source_passage
from .step_selection import SOURCE_OPERATIONS, _call, _selection

QuestionStatus = Literal["open", "investigating", "resolved", "withdrawn"]
Reference = Annotated[str, Field(pattern=r"^sha256:[0-9a-f]{64}$")]
Text = Annotated[str, Field(min_length=1, max_length=2000)]


class ScientificQuestion(BaseModel):
    """A scoped agent judgment; resolution does not establish chemical validity."""

    model_config = ConfigDict(extra="forbid", frozen=True, str_strip_whitespace=True)
    schema_version: Literal["scientific_question.v1"] = "scientific_question.v1"
    question_id: str
    source_ref: Reference
    step_id: Text | None = None
    question: Text
    resolution_criterion: Text
    next_check: Text
    status: QuestionStatus = "open"
    rationale: Text | None = None
    evidence_refs: tuple[Reference, ...] = Field(default=(), max_length=20)
    contradicting_refs: tuple[Reference, ...] = Field(default=(), max_length=20)
    resolution_refs: tuple[Reference, ...] = Field(default=(), max_length=20)
    previous_ref: Reference | None = None
    origin: Literal["agent_authored"] = "agent_authored"
    review_status: Literal["unreviewed"] = "unreviewed"

    @model_validator(mode="after")
    def validate_resolution(self) -> ScientificQuestion:
        """Require an explanation for revisions and explicit evidence for resolution."""
        if self.previous_ref is not None and not self.rationale:
            raise ValueError("Question updates require a rationale")
        if self.status == "resolved" and (not self.resolution_refs or not self.rationale):
            raise ValueError("A resolved question requires resolution_refs and a rationale")
        if self.status != "resolved" and self.resolution_refs:
            raise ValueError("resolution_refs are only valid for resolved questions")
        return self


def _heads(store: InvestigationStore) -> dict[str, tuple[str, ScientificQuestion]]:
    heads: dict[str, tuple[str, ScientificQuestion]] = {}
    for event in store.events():
        if event.kind != "scientific_question":
            continue
        question = ScientificQuestion.model_validate(store.read_artifact(event.artifact_ref))
        previous = heads.get(question.question_id)
        if question.previous_ref != (previous[0] if previous else None):
            raise ValueError("Question history is not a single revision chain")
        heads[question.question_id] = (event.artifact_ref, question)
    return heads


def _references(store: InvestigationStore, refs: tuple[str, ...], *, resolution: bool = False) -> None:
    events = {event.artifact_ref: event for event in store.events()}
    for ref in refs:
        event = events.get(ref)
        if event is None:
            raise ValueError("Question evidence must identify a recorded event")
        payload = store.read_artifact(ref)
        if resolution:
            if event.kind in {"literature_source", "literature_excerpt"}:
                load_source_passage(store, ref)
            elif event.kind != "call" or payload.get("execution_status") != "completed":
                raise ValueError("Resolution needs a completed scientific call or captured source passage")
        elif event.kind not in {"call", "literature_source", "literature_excerpt"}:
            raise ValueError("Use scientific calls or source records, not agent notes, as question evidence")


def record_question(
    store: InvestigationStore, *, source_ref: str, question: str,
    resolution_criterion: str, next_check: str, step_id: str | None = None,
) -> InvestigationEvent:
    """Create one question about an exact saved call or selected assessed route step."""
    with store._writer():
        _references(store, (source_ref,))
        if step_id is not None:
            payload = _call(store, source_ref, SOURCE_OPERATIONS - {"disconnect_target"})
            _selection(payload, step_id, None)
        value = ScientificQuestion(
            question_id="pending", source_ref=source_ref, step_id=step_id,
            question=question, resolution_criterion=resolution_criterion, next_check=next_check,
        )
        identity = {name: getattr(value, name) for name in (
            "schema_version", "source_ref", "step_id", "question", "resolution_criterion",
        )}
        question_id = "q_" + hashlib.sha256(canonical_bytes(identity)).hexdigest()
        if question_id in _heads(store):
            raise ValueError("Question already exists; update its current revision")
        value = value.model_copy(update={"question_id": question_id})
        return store._append_locked("scientific_question", value.model_dump(mode="json"),
                                    evidence_refs=(source_ref,))


def update_question(
    store: InvestigationStore, reference: str, *, status: QuestionStatus, rationale: str,
    next_check: str | None = None, evidence_refs: tuple[str, ...] = (),
    contradicting_refs: tuple[str, ...] = (), resolution_refs: tuple[str, ...] = (),
) -> InvestigationEvent:
    """Append to the current revision, retaining contradictions and prior resolutions."""
    with store._writer():
        head = next((item for item in _heads(store).values() if item[0] == reference), None)
        if head is None:
            raise ValueError("Update requires the current question revision_ref; refresh list_questions")
        previous = head[1]
        value = ScientificQuestion.model_validate({
            **previous.model_dump(), "previous_ref": reference, "status": status,
            "rationale": rationale, "next_check": next_check if next_check is not None else previous.next_check,
            "evidence_refs": tuple(dict.fromkeys((*previous.evidence_refs, *evidence_refs))),
            "contradicting_refs": tuple(dict.fromkeys((*previous.contradicting_refs, *contradicting_refs))),
            "resolution_refs": resolution_refs,
        })
        _references(store, (*value.evidence_refs, *value.contradicting_refs))
        if status == "resolved":
            _references(store, value.resolution_refs, resolution=True)
        refs = tuple(dict.fromkeys((value.source_ref, reference, *value.evidence_refs,
                                   *value.contradicting_refs, *value.resolution_refs)))
        return store._append_locked("scientific_question", value.model_dump(mode="json"), evidence_refs=refs)


def list_questions(
    store: InvestigationStore, *, source_ref: str | None = None, step_id: str | None = None,
    status: QuestionStatus | None = None, offset: int = 0, limit: int = 20,
) -> dict[str, Any]:
    """Page current question revisions; exact scope filters never transfer a resolution."""
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 20:
        raise ValueError("Use nonnegative offset and limit between 1 and 20")
    if status is not None and status not in {"open", "investigating", "resolved", "withdrawn"}:
        raise ValueError("Unknown question status")
    if step_id is not None and source_ref is None:
        raise ValueError("step_id requires source_ref")
    if source_ref is not None:
        _references(store, (source_ref,))
    items = [{"revision_ref": ref, **question.model_dump(mode="json")}
             for ref, question in _heads(store).values()
             if (source_ref is None or (question.source_ref == source_ref and question.step_id == step_id))
             and (status is None or question.status == status)]
    return {"schema_version": "scientific_question_list.v1", "items": items[offset:offset + limit],
            "total": len(items), "offset": offset,
            "next_offset": offset + limit if offset + limit < len(items) else None,
            "resolution_authority": "agent_assessed_unreviewed_not_chemical_validation"}
