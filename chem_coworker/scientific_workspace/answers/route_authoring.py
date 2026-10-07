"""Compile compact route input into the canonical answer view without running science."""

from __future__ import annotations

from typing import Any, Literal, Mapping

from pydantic import Field

from reactive_taxonomy.reaction_operators import canonical_molecule_collection

from ..adapters.literature import load_source_passage
from ..core.store import InvestigationStore
from .answer_contracts import AnswerObject, ScientificAnswer
from .answer_finalization import _complete_empty_fields


class RouteSupport(AnswerObject):
    """An inspected saved source; URLs and titles come from the evidence store."""

    ref: str = Field(pattern=r"^sha256:[0-9a-f]{64}$")
    locator: str = Field(min_length=1, max_length=1000)
    note: str = Field(default="", max_length=2000)


class RouteStepInput(AnswerObject):
    """One proposed or reported reaction, with optional conditions and support."""

    reaction_smiles: str = Field(min_length=1, max_length=12000)
    title: str = Field(default="", max_length=300)
    basis: Literal["proposed", "reported"] = "proposed"
    conditions: str = Field(default="", max_length=4000)
    conditions_basis: Literal["proposed", "reported", "unknown"] = "proposed"
    reagents: list[str] = Field(default_factory=list, max_length=30)
    support: list[RouteSupport] = Field(default_factory=list, max_length=30)
    limitations: list[str] = Field(default_factory=list, max_length=30)
    after_steps: list[int] | None = Field(default=None, max_length=20)


class RouteInput(AnswerObject):
    """Steps in synthetic order; one-based dependency overrides stay route-local."""

    title: str = Field(default="", max_length=300)
    steps: list[RouteStepInput] = Field(min_length=1, max_length=30)
    limitations: list[str] = Field(default_factory=list, max_length=30)


class RouteAnswerInput(AnswerObject):
    """Versioned authoring input; stored answers still use scientific_answer.v2."""

    schema_version: Literal["route_answer_input.v1"] = "route_answer_input.v1"
    target_smiles: str = Field(min_length=1, max_length=4000)
    routes: list[RouteInput] = Field(min_length=1, max_length=10)
    summary: str = Field(default="Proposed routes, conditions and supporting literature.", min_length=1)
    limitations: list[str] = Field(default_factory=list, max_length=30)


def _components(smiles: str) -> list[str]:
    """Delegate graph identity, stereo and map normalization to the taxonomy."""
    if not smiles.strip() or any(not part.strip() for part in smiles.split(".")):
        raise ValueError("Supply nonempty molecular components on both reaction sides")
    canonical = canonical_molecule_collection(smiles)
    if not canonical:
        raise ValueError(f"Invalid route SMILES: {smiles[:120]}")
    return canonical.split(".")


def route_answer(store: InvestigationStore, proposal: Mapping[str, Any]) -> dict[str, Any]:
    """Assemble IDs, dependencies and citations; never infer conditions or feasibility.

    Only the supplied sources are attached. Existing inspection references retain
    the canonical finalizer's exact-structure checks. Missing required inspections
    still fail publication. No retrieval, assessment, source reconstruction or
    self-review is performed. Canonicalization identifies display molecules only.
    """
    request = RouteAnswerInput.model_validate(dict(proposal))
    draft = _complete_empty_fields({"answer_markdown": request.summary,
                                    "uncertainties": request.limitations})
    molecule_ids: dict[str, str] = {}
    source_ids: dict[tuple[str, str], str] = {}
    kinds = {event.artifact_ref: event.kind for event in store.events()}

    def molecule(smiles: str, *, target: bool = False) -> str:
        if smiles not in molecule_ids:
            identity = f"m{len(molecule_ids) + 1}"
            molecule_ids[smiles] = identity
            draft["molecules"].append({"id": identity, "name": "Target" if target else identity,
                                       "smiles": smiles, "basis": "input" if target else "proposed"})
        return molecule_ids[smiles]

    def source(support: RouteSupport) -> tuple[str, str | None]:
        key = (support.ref, support.locator)
        operation = None
        kind = kinds.get(support.ref)
        if kind in {"literature_source", "literature_excerpt"}:
            saved, _, _ = load_source_passage(store, support.ref)
            citation = {"kind": "external_source", "url": saved.get("final_url") or saved["source_url"],
                        "title": saved.get("title") or support.locator}
        elif kind == "call":
            saved = store.read_artifact(support.ref)
            if saved.get("execution_status") != "completed":
                raise ValueError("Route support requires a completed recorded call")
            operation = saved["operation"]
            citation = {"kind": "local_artifact", "url": None, "title": operation.replace("_", " ")}
        else:
            raise ValueError("Route support requires captured literature or a recorded scientific call")
        if key not in source_ids:
            identity = f"source{len(source_ids) + 1}"
            source_ids[key] = identity
            draft["sources"].append({"id": identity, **citation, "artifact_ref": support.ref,
                                     "locator": support.locator})
        return source_ids[key], operation

    draft["target_molecule_ids"] = [molecule(part, target=True) for part in _components(request.target_smiles)]
    for route_number, route in enumerate(request.routes, 1):
        route_id = f"r{route_number}"
        step_ids = []
        producers: dict[str, list[int]] = {}
        for step_number, step in enumerate(route.steps, 1):
            parts = step.reaction_smiles.split(">")
            if len(parts) != 3 or parts[1]:
                raise ValueError("Use reactants>>products; put reagents and conditions in their separate fields")
            left, right = _components(parts[0]), _components(parts[2])
            reactants = [molecule(part) for part in left]
            products = [molecule(part) for part in right]
            if step.after_steps is None:
                if any(len(producers.get(part, [])) > 1 for part in left):
                    raise ValueError("Ambiguous intermediate producer; specify one-based after_steps for this route")
                parents = sorted({number for part in left for number in producers.get(part, [])})
            else:
                parents = step.after_steps
                if len(set(parents)) != len(parents) or any(number < 1 or number >= step_number for number in parents):
                    raise ValueError("after_steps must name distinct earlier one-based steps in this route")
            identity = f"{route_id}_s{step_number}"
            citations, notes = [], list(step.limitations)
            links: dict[str, list[str]] = {key: [] for key in (
                "precedent_refs", "condition_precedent_refs", "recipe_assessment_refs",
            )}
            link_fields = {"inspect_step_precedents": "precedent_refs",
                           "inspect_condition_precedents": "condition_precedent_refs",
                           "assess_proposed_recipe": "recipe_assessment_refs"}
            for support in step.support:
                citation_id, operation = source(support)
                if citation_id not in citations:
                    citations.append(citation_id)
                if support.note:
                    notes.append(f"{support.locator}: {support.note}")
                field = link_fields.get(operation)
                if field and support.ref not in links[field]:
                    links[field].append(support.ref)
            condition_fields = {"basis": step.conditions_basis, "source_ids": citations}
            draft["steps"].append({
                "id": identity, "title": step.title or f"Step {step_number}",
                "basis": step.basis, "source_ids": citations, "limitations": notes,
                "reactant_ids": reactants, "product_ids": products,
                "after_step_ids": [f"{route_id}_s{number}" for number in parents],
                "conditions": [{"text": step.conditions, **condition_fields}] if step.conditions else [],
                "reagents": [{"text": text, **condition_fields} for text in step.reagents],
                **links,
            })
            step_ids.append(identity)
            for part in set(right):
                producers.setdefault(part, []).append(step_number)
        draft["routes"].append({"id": route_id, "title": route.title or f"Route {route_number}",
                                 "step_ids": step_ids, "limitations": route.limitations})
    return ScientificAnswer.model_validate(_complete_empty_fields(draft)).model_dump()
