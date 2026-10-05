"""Small immutable authoring helpers; no chemistry inference or automatic tool calls."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping

from ..adapters.literature import load_source_passage
from ..adapters.literature_reactions import prepared_literature_reaction
from ..core.store import InvestigationStore
from .answer_contracts import ScientificAnswer, validate_answer_evidence
from .answer_finalization import _complete_empty_fields


def attach_literature_reaction(
    store: InvestigationStore, draft: Mapping[str, Any], step_id: str, preparation_ref: str,
) -> dict[str, Any]:
    """Attach an unchanged preparation and its exact required AnswerSource together."""
    value = deepcopy(dict(draft))
    steps = [step for step in value.get("steps", []) if step.get("id") == step_id]
    if len(steps) != 1:
        raise ValueError("Select one actual answer step_id")
    reaction = prepared_literature_reaction(store, preparation_ref)
    record = store.read_artifact(preparation_ref)["result"]
    source, _, _ = load_source_passage(store, record["source_ref"])
    citation = {"id": reaction["source_id"], "kind": "external_source", "title": reaction["title"],
                "artifact_ref": record["source_ref"], "url": source.get("final_url") or source["source_url"],
                "locator": reaction["locator"]}
    existing = [item for item in value.setdefault("sources", []) if item["id"] == citation["id"]]
    if existing and (len(existing) != 1 or existing[0]["artifact_ref"] != citation["artifact_ref"]
                     or existing[0].get("url") != citation["url"]):
        raise ValueError("Existing source_id conflicts with preparation; use a distinct source_id during preparation")
    if not existing:
        value["sources"].append(citation)
    reactions = steps[0].setdefault("literature_reactions", [])
    if not any(item.get("preparation_ref") == preparation_ref for item in reactions):
        reactions.append(reaction)
    return value


def attach_recipe_check(
    store: InvestigationStore, draft: Mapping[str, Any], step_id: str, reference: str,
) -> dict[str, Any]:
    """Attach only a saved check for the exact answer-step structures and stereo."""
    from ..adapters.investigation_checks import load_proposed_recipe_check

    value = deepcopy(dict(draft))
    molecules = {item["id"]: item["smiles"] for item in value.get("molecules", [])}
    steps = [step for step in value.get("steps", []) if step.get("id") == step_id]
    if len(steps) != 1:
        raise ValueError("Select one actual answer step_id")
    step = steps[0]
    load_proposed_recipe_check(store, reference, ".".join(molecules[key] for key in step["reactant_ids"]),
                              ".".join(molecules[key] for key in step["product_ids"]))
    refs = step.setdefault("recipe_assessment_refs", [])
    if reference not in refs:
        refs.append(reference)
    return value


def answer_preflight(
    store: InvestigationStore, draft: Mapping[str, Any], offset: int = 0, limit: int = 10,
) -> dict[str, Any]:
    """Report validation errors and missing recorded checks without running science."""
    warnings = []
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 10:
        raise ValueError("offset must be nonnegative; limit must be 1..10")
    try:
        answer = ScientificAnswer.model_validate(_complete_empty_fields(draft))
        refs = validate_answer_evidence(answer, store)
    except (ValueError, TypeError, OSError, KeyError) as exc:
        return {"schema_version": "answer_preflight.v1", "valid": False,
                "error": str(exc)[:3000], "warnings": [], "science_rerun": False}
    for step in answer.steps:
        for reaction in step.literature_reactions:
            if reaction.preparation_ref:
                preparation = store.read_artifact(reaction.preparation_ref)["result"]
                conflicts = [item for item in preparation.get("quantity_checks", []) if item["status"] == "conflicting"]
                if conflicts:
                    warnings.append({"step_id": step.id, "gap": "source_quantity_conflict",
                                     "artifact_ref": reaction.preparation_ref, "conflict_count": len(conflicts),
                                     "action": "Inspect both reported quantities; disclose the discrepancy and do not silently derive a recipe from one value."})
        if step.conditions and not step.recipe_assessment_refs:
            warnings.append({"step_id": step.id, "gap": "actual_recipe_check_not_attached",
                             "action": "Normalize the actual recipe with assess_proposed_recipe and attach_recipe_check; disclose unresolved coverage."})
        for ref in step.recipe_assessment_refs:
            saved = store.read_artifact(ref)["result"]
            assessment = saved["compatibility"]
            completeness = assessment.get("reaction_completeness") or {}
            if completeness.get("product_element_excess"):
                warnings.append({"step_id": step.id, "gap": "reaction_inputs_incomplete",
                                 "artifact_ref": ref, "product_element_excess": completeness["product_element_excess"],
                                 "action": "Inspect evidence-backed atom contributors and stoichiometric multiplicity in the reaction graph. "
                                 "Recipe quantities do not supply graph atoms. Do not invent donors or mapping."})
            if (assessment.get("status") != "no_known_conflict"
                    or assessment.get("coverage", {}).get("capability_status") != "supported"
                    or saved.get("process_coverage") == "stages_recorded_not_evaluated"):
                warnings.append({"step_id": step.id, "gap": "recipe_coverage_or_conflicts_unresolved",
                                 "artifact_ref": ref, "status": assessment.get("status"),
                                 "action": "Read coverage and hard_conflicts; retain uncertainty. A valid citation does not establish recipe support."})
    route_refs = []
    input_reviews = {}
    indirect_routes = set()
    for event in store.events():
        if event.kind != "call":
            continue
        payload = store.read_artifact(event.artifact_ref)
        if payload.get("execution_status") != "completed":
            continue
        if payload.get("operation") == "inspect_route_inputs":
            input_reviews[payload["result"]["source_ref"]] = (event.artifact_ref, payload["result"])
        if event.artifact_ref in refs and payload.get("operation") in {"inspect_step_precedents", "inspect_route_step"}:
            indirect_routes.add(payload["result"]["source_ref"])
        if event.artifact_ref in refs and payload.get("operation") in {"assess_route_proposal", "revise_route_branch"}:
            route_refs.append(event.artifact_ref)
    for ref in indirect_routes:
        record = store.read_artifact(ref)
        if record.get("operation") in {"assess_route_proposal", "revise_route_branch"}:
            route_refs.append(ref)
    for ref in dict.fromkeys(route_refs):
        if ref not in input_reviews:
            warnings.append({"gap": "route_inputs_not_inspected", "source_ref": ref,
                             "action": "Use inspect_route_inputs and search captured sources for upstream preparations; disclose advanced-input scope."})
        else:
            review_ref, review = input_reviews[ref]
            unresolved = [item["smiles"] for item in review["leaves"]
                          if item["starting_material_assessment"].get("status") in {"unresolved", "assumed_terminal"}
                          or item["source_search_status"] != "searched"]
            if unresolved:
                warnings.append({"gap": "route_inputs_remain_assumed_or_unsearched", "artifact_ref": review_ref,
                                 "leaf_count": len(unresolved),
                                 "action": "Read leaf status and captured upstream leads; disclose supply assumptions and the route's starting scope."})
    if offset > len(warnings):
        raise ValueError("offset exceeds the preflight warning count")
    return {"schema_version": "answer_preflight.v1", "valid": True, "evidence_count": len(refs),
            "warnings": warnings[offset:offset + limit], "warning_count": len(warnings),
            "next_offset": offset + limit if offset + limit < len(warnings) else None, "science_rerun": False,
            "limitations": ["Validation checks citations and structure association, not scientific truth or experimental feasibility."]}
