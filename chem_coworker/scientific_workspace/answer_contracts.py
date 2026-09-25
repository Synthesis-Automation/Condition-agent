"""Agent-authored scientific answer views, separate from authoritative domain records."""

from __future__ import annotations

import re
from typing import Any, Literal
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .store import InvestigationStore


Basis = Literal["reported", "computed", "proposed", "unknown", "input"]


class AnswerObject(BaseModel):
    """Closed JSON objects supplied by an agent, never admitted chemistry records."""

    model_config = ConfigDict(extra="forbid", strict=True)


class AnswerSource(AnswerObject):
    """A precise local result or captured external-source excerpt for inspection."""

    id: str = Field(pattern=r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")
    kind: Literal["local_artifact", "external_source"]
    title: str = Field(min_length=1, max_length=500)
    artifact_ref: str = Field(pattern=r"^sha256:[0-9a-f]{64}$")
    url: str | None
    locator: str = Field(min_length=1, max_length=1000)

    @model_validator(mode="after")
    def check_location(self) -> "AnswerSource":
        """External references require an inspectable web location and saved excerpt."""
        if self.kind == "external_source":
            if not self.url or urlsplit(self.url).scheme not in {"http", "https"} or not urlsplit(self.url).netloc:
                raise ValueError("External sources require an HTTP(S) URL")
        elif self.url is not None:
            raise ValueError("Local artifact sources must use url=null")
        return self


class AttributedObject(AnswerObject):
    """Explicit epistemic status and citations for a statement or scientific object."""

    basis: Basis
    source_ids: list[str] = Field(max_length=30)
    limitations: list[str] = Field(max_length=30)

    @model_validator(mode="after")
    def require_support(self) -> "AttributedObject":
        """Reported and computed labels require citations, not just confident wording."""
        if self.basis in {"reported", "computed"} and not self.source_ids:
            raise ValueError("Reported and computed objects require source_ids")
        return self


class AnswerClaim(AttributedObject):
    """A condition, yield, or standalone statement with its own attribution."""

    text: str = Field(min_length=1, max_length=4000)


class AnswerMolecule(AttributedObject):
    """A named structure proposed for display; SMILES validity is not implied."""

    id: str = Field(pattern=r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")
    name: str = Field(min_length=1, max_length=300)
    smiles: str = Field(min_length=1, max_length=4000)


class AnswerStep(AttributedObject):
    """One declared transformation with separately attributed conditions and yield."""

    id: str = Field(pattern=r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")
    title: str = Field(min_length=1, max_length=300)
    reactant_ids: list[str] = Field(min_length=1, max_length=20)
    product_ids: list[str] = Field(min_length=1, max_length=20)
    after_step_ids: list[str] = Field(max_length=20)
    conditions: list[AnswerClaim] = Field(max_length=30)
    yield_info: AnswerClaim | None


class AnswerRoute(AnswerObject):
    """An explicitly ordered route view; feasibility remains a separate assessment."""

    id: str = Field(pattern=r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")
    title: str = Field(min_length=1, max_length=300)
    step_ids: list[str] = Field(min_length=1, max_length=30)
    limitations: list[str] = Field(max_length=30)


class ScientificAnswer(AnswerObject):
    """Versioned, unreviewed answer contract for new turns; legacy saved text stays readable."""

    schema_version: Literal["scientific_answer.v2"]
    answer_markdown: str = Field(min_length=1)
    evidence_refs: list[str]
    uncertainties: list[str]
    needs_user_input: bool
    sources: list[AnswerSource] = Field(max_length=40)
    molecules: list[AnswerMolecule] = Field(max_length=60)
    target_molecule_ids: list[str] = Field(max_length=10)
    steps: list[AnswerStep] = Field(max_length=30)
    routes: list[AnswerRoute] = Field(max_length=10)
    claims: list[AnswerClaim] = Field(max_length=40)

    def attributed_objects(self) -> list[AttributedObject]:
        """Visit every object whose status/citations must be checked."""
        objects: list[AttributedObject] = [*self.molecules, *self.steps, *self.claims]
        for step in self.steps:
            objects.extend(step.conditions)
            if step.yield_info is not None:
                objects.append(step.yield_info)
        return objects

    @model_validator(mode="after")
    def check_references(self) -> "ScientificAnswer":
        """Check view identity and dependency integrity without inventing chemistry."""
        for group in (self.sources, self.molecules, self.steps, self.routes):
            if len({item.id for item in group}) != len(group):
                raise ValueError("Scientific object IDs must be unique within each collection")
        sources = {item.id: item for item in self.sources}
        molecules = {item.id for item in self.molecules}
        steps = {item.id: item for item in self.steps}
        for item in self.attributed_objects():
            if not set(item.source_ids) <= sources.keys():
                raise ValueError("Unknown source_id in scientific answer")
            if item.basis == "computed" and not any(sources[key].kind == "local_artifact" for key in item.source_ids):
                raise ValueError("Computed objects require a local computation source")
        if not set(self.target_molecule_ids) <= molecules:
            raise ValueError("Unknown target molecule ID")
        visited: set[str] = set()
        for step in self.steps:
            if not set(step.reactant_ids + step.product_ids) <= molecules:
                raise ValueError("Unknown molecule ID in reaction step")
            if not set(step.after_step_ids) <= visited:
                raise ValueError("Steps must be acyclic and ordered after their dependencies")
            for parent in step.after_step_ids:
                if not set(steps[parent].product_ids).intersection(step.reactant_ids):
                    raise ValueError("Dependent steps must share an explicit intermediate molecule ID")
            visited.add(step.id)
        for route in self.routes:
            included: set[str] = set()
            for step_id in route.step_ids:
                if step_id not in steps or step_id in included:
                    raise ValueError("Unknown or duplicate route step ID")
                if not set(steps[step_id].after_step_ids) <= included:
                    raise ValueError("Route must include each step's preceding dependencies")
                included.add(step_id)
        return self


ANSWER_SCHEMA = ScientificAnswer.model_json_schema()


def validate_answer_evidence(answer: ScientificAnswer, store: InvestigationStore) -> list[str]:
    """Verify recorded citations and their types; do not claim semantic fact-checking."""
    kinds = {event.artifact_ref: event.kind for event in store.events()}
    cited = set(answer.evidence_refs)
    cited.update(re.findall(r"sha256:[0-9a-f]{64}", answer.model_dump_json()))
    for reference in cited:
        store.read_artifact(reference)
        if kinds.get(reference) not in {"call", "derived_file", "replay", "custom_execution"}:
            raise ValueError("Answer must cite scientific evidence, not agent assertions")
    sources = {source.id: source for source in answer.sources}
    for source in answer.sources:
        if source.kind == "external_source" and kinds[source.artifact_ref] != "derived_file":
            raise ValueError("External sources require a saved source excerpt attachment")
    for item in answer.attributed_objects():
        if item.basis != "computed":
            continue
        supports: list[tuple[str, Any]] = [
            (kinds[sources[key].artifact_ref], store.read_artifact(sources[key].artifact_ref))
            for key in item.source_ids if sources[key].kind == "local_artifact"
        ]
        if not any(isinstance(value, dict) and (
            kind == "call" and value.get("execution_status") == "completed" and "operation" in value
            or kind == "replay" and "source_ref" in value and "matches" in value
            or kind == "custom_execution" and value.get("execution_status") == "completed"
            and value.get("returncode") == 0 and "script_sha256" in value and "input_sha256" in value
        ) for kind, value in supports):
            raise ValueError("Computed objects require a completed recorded computation or replay")
    return sorted(cited)
