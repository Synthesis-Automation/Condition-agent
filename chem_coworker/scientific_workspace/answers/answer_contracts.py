"""Agent-authored scientific answer views, separate from authoritative domain records."""

from __future__ import annotations

import re
from typing import Any, Literal
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..core.store import InvestigationStore

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


class LiteratureStructure(AnswerObject):
    """An explicitly supplied source structure or declared reconstruction."""

    name: str = Field(min_length=1, max_length=300)
    smiles: str = Field(min_length=1, max_length=4000)


class LiteratureReaction(AnswerObject):
    """Source-linked drawing input, never an admitted or indexed precedent."""

    schema_version: Literal["literature_reaction.v1"] = "literature_reaction.v1"
    title: str = Field(min_length=1, max_length=300)
    source_id: str
    locator: str = Field(min_length=1, max_length=1000)
    structure_origin: Literal["source_explicit", "reconstructed_from_description"]
    structure_evidence: str = Field(min_length=1, max_length=2000)
    reactants: list[LiteratureStructure] = Field(min_length=1, max_length=20)
    products: list[LiteratureStructure] = Field(min_length=1, max_length=20)
    conditions: list[AnswerClaim] = Field(default_factory=list, max_length=30)
    yield_info: AnswerClaim | None = None
    limitations: list[str] = Field(default_factory=list, max_length=30)


class AnswerStep(AttributedObject):
    """A transformation with attributed reagent labels, procedures and yield."""

    id: str = Field(pattern=r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")
    title: str = Field(min_length=1, max_length=300)
    reactant_ids: list[str] = Field(min_length=1, max_length=20)
    product_ids: list[str] = Field(min_length=1, max_length=20)
    after_step_ids: list[str] = Field(max_length=20)
    conditions: list[AnswerClaim] = Field(max_length=30)
    reagents: list[AnswerClaim] = Field(default_factory=list, max_length=30)
    yield_info: AnswerClaim | None
    precedent_refs: list[str] = Field(default_factory=list, max_length=5)
    condition_precedent_refs: list[str] = Field(default_factory=list, max_length=5)
    rationale: AnswerClaim | None = None
    literature_reactions: list[LiteratureReaction] = Field(default_factory=list, max_length=5)


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
            objects.extend(step.reagents)
            if step.rationale is not None:
                objects.append(step.rationale)
            if step.yield_info is not None:
                objects.append(step.yield_info)
            for reaction in step.literature_reactions:
                objects.extend(reaction.conditions)
                if reaction.yield_info is not None:
                    objects.append(reaction.yield_info)
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
            for reaction in step.literature_reactions:
                source = sources.get(reaction.source_id)
                if source is None or source.kind != "external_source":
                    raise ValueError("Literature drawings require a captured external source")
                for claim in [*reaction.conditions, *([reaction.yield_info] if reaction.yield_info else [])]:
                    if claim.basis != "reported" or reaction.source_id not in claim.source_ids:
                        raise ValueError("Literature conditions and yields must be reported by the drawing's source")
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
# Runtime schemas declare every field explicitly (including nullable rationale).
# Older saved v2 answers still load through model defaults without invented support.
ANSWER_SCHEMA["$defs"]["AnswerStep"]["required"].extend([
    "precedent_refs", "condition_precedent_refs", "rationale", "reagents", "literature_reactions",
])


def validate_answer_evidence(answer: ScientificAnswer, store: InvestigationStore) -> list[str]:
    """Verify recorded citations and their types; do not claim semantic fact-checking."""
    kinds = {event.artifact_ref: event.kind for event in store.events()}
    cited = set(answer.evidence_refs)
    cited.update(re.findall(r"sha256:[0-9a-f]{64}", answer.model_dump_json()))
    for reference in cited:
        store.read_artifact(reference)
        if kinds.get(reference) not in {
            "call", "derived_file", "replay", "custom_execution", "literature_source", "literature_excerpt",
        }:
            raise ValueError("Answer must cite scientific evidence, not agent assertions")
    from .condition_precedents import load_condition_precedent_evidence
    from .step_precedents import (
        load_step_precedent_evidence,
        require_available_step_precedents,
    )

    molecules = {item.id: item.smiles for item in answer.molecules}
    for step in answer.steps:
        for reference in step.precedent_refs:
            load_step_precedent_evidence(
                store, reference, ".".join(molecules[key] for key in step.reactant_ids),
                ".".join(molecules[key] for key in step.product_ids),
            )
        for reference in step.condition_precedent_refs:
            load_condition_precedent_evidence(
                store, reference, ".".join(molecules[key] for key in step.reactant_ids),
                ".".join(molecules[key] for key in step.product_ids),
            )
    require_available_step_precedents(store, answer.model_dump())
    sources = {source.id: source for source in answer.sources}
    for source in answer.sources:
        if source.kind != "external_source":
            continue
        kind = kinds[source.artifact_ref]
        if kind not in {"derived_file", "literature_source", "literature_excerpt"}:
            raise ValueError("External sources require a saved source excerpt attachment or literature snapshot")
        if kind in {"literature_source", "literature_excerpt"}:
            value = store.read_artifact(source.artifact_ref)
            text = value.get("text") if kind == "literature_excerpt" else value.get("extraction", {}).get("text")
            if not isinstance(text, str) or not text.strip():
                raise ValueError("External sources require captured text; failed retrieval is not a source passage")
            def location(url: str) -> str:
                parsed = urlsplit(url)
                return parsed._replace(fragment="", path=parsed.path or "/").geturl()

            recorded_urls = {location(url) for url in (value.get("source_url"), value.get("final_url")) if url}
            if location(source.url or "") not in recorded_urls:
                raise ValueError("External source URL does not match its captured source")
    for step in answer.steps:
        for reaction in step.literature_reactions:
            source = sources[reaction.source_id]
            kind = kinds[source.artifact_ref]
            if kind not in {"literature_source", "literature_excerpt"}:
                raise ValueError("Literature drawings require recorded captured text, not an arbitrary attachment")
            value = store.read_artifact(source.artifact_ref)
            text = value.get("text") if kind == "literature_excerpt" else value.get("extraction", {}).get("text", "")
            if " ".join(reaction.structure_evidence.split()) not in " ".join(text.split()):
                raise ValueError("Literature drawing structure_evidence must quote its captured source")
            if reaction.structure_origin == "source_explicit":
                if any(item.smiles not in text for item in [*reaction.reactants, *reaction.products]):
                    raise ValueError("Source-explicit SMILES must appear in captured text; otherwise label a reconstruction")
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
