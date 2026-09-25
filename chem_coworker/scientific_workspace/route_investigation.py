"""Recorded application composition for external route investigation and revision."""

from __future__ import annotations

from dataclasses import fields
from typing import Any, TYPE_CHECKING

from core_retrosynthesis.external_proposal_assessment import (
    ExternalRetrosynthesisProposal, assess_external_retrosynthesis_proposal,
)
from core_retrosynthesis.external_route_admission import (
    ExternalRouteProposal, ExternalRouteStepProposal, assess_external_route_proposal,
    external_route_proposal_from_tree,
)
from core_retrosynthesis.route_proposal_revision import (
    assess_declared_route_materials, external_route_step_neighbors, revise_external_route_proposal,
)

if TYPE_CHECKING:
    from .operations import ScientificOperations


ROUTE_RECORD_OPERATIONS = frozenset({"assess_route_proposal", "prepare_route_proposal", "revise_route_branch"})


def _strings(values: list[str] | None, name: str) -> list[str]:
    if values is None:
        return []
    if not isinstance(values, list) or not all(isinstance(item, str) and item.strip() for item in values):
        raise ValueError(f"{name} must be a list of nonempty strings")
    return list(values)


def _evidence(operations: ScientificOperations, references: list[str] | None) -> list[str]:
    refs = _strings(references, "evidence_refs")
    kinds = {event.artifact_ref: event.kind for event in operations.store.events()}
    for ref in refs:
        operations.store.read_artifact(ref)
        if kinds.get(ref) not in {"call", "derived_file", "replay", "custom_execution"}:
            raise ValueError("Route evidence must reference recorded scientific evidence")
    return list(dict.fromkeys(refs))


def _known_fields(value: dict[str, Any], contract: type, extra: tuple[str, ...] = ()) -> None:
    if not isinstance(value, dict) or set(value) - {field.name for field in fields(contract)} - set(extra):
        raise ValueError(f"Unknown fields or invalid object for {contract.__name__}")


def _step(value: dict[str, Any], *, named: bool = False) -> Any:
    _known_fields(value, ExternalRetrosynthesisProposal, ("external_step_id",) if named else ())
    for key in ("target_smiles", "precursor_smiles", *(("external_step_id",) if named else ())):
        if not isinstance(value.get(key), str) or not value[key].strip():
            raise ValueError(f"{key} must be a nonempty string")
    return ExternalRouteStepProposal.from_dict(value) if named else ExternalRetrosynthesisProposal.from_dict(value)


def _proposal(value: dict[str, Any]) -> ExternalRouteProposal:
    _known_fields(value, ExternalRouteProposal)
    if not isinstance(value.get("target_smiles"), str) or not value["target_smiles"].strip():
        raise ValueError("target_smiles must be a nonempty string")
    if not isinstance(value.get("steps"), list) or len(value["steps"]) > 40:
        raise ValueError("steps must be a list of at most 40 physical steps")
    for step in value["steps"]:
        _step(step, named=True)
    return ExternalRouteProposal.from_dict(value)


def _assessment_arguments(operations: ScientificOperations, include_conditions: bool, include_forward: bool) -> dict[str, Any]:
    if type(include_conditions) is not bool or type(include_forward) is not bool:
        raise ValueError("Assessment options must be booleans")
    library = operations._external_route_library()
    arguments: dict[str, Any] = {"operator_library": library}
    if include_forward:
        from core_retrosynthesis.forward_assessment import build_forward_library_from_generic

        arguments["forward_library"] = build_forward_library_from_generic(library)
    if include_conditions:
        from core_retrosynthesis import recommend_retrosynthesis_conditions

        recommender = operations._conditions()
        arguments["condition_evaluator"] = lambda reaction: recommend_retrosynthesis_conditions(
            reaction, recommender, condition_top_k=3,
        )
    return arguments


def _recipe(operations: ScientificOperations, proposal: ExternalRetrosynthesisProposal) -> Any:
    # Assess a supplied recipe separately from retrieved analogue conditions.
    if proposal.proposed_conditions is None:
        return None
    return operations.assess_recipe(
        f"{proposal.precursor_smiles}>>{proposal.target_smiles}", dict(proposal.proposed_conditions),
    )


def assess_step(
    operations: ScientificOperations, proposal: dict[str, Any], include_conditions: bool,
    include_forward: bool, evidence_refs: list[str] | None,
) -> dict[str, Any]:
    """Compose canonical single-step assessment and supplied-recipe checks."""
    refs = _evidence(operations, evidence_refs)
    step = _step(proposal)
    assessment = assess_external_retrosynthesis_proposal(
        step, **_assessment_arguments(operations, include_conditions, include_forward),
    )
    return {
        "schema_version": "route_step_investigation.v1", "origin": "agent_proposal",
        "review_status": "unreviewed", "proposal": step.to_dict(),
        "assessment": assessment.to_dict(), "proposed_recipe_assessment": _recipe(operations, step),
        "evidence_refs": refs, "experimental_feasibility": "not_established",
    }


def assess_route(
    operations: ScientificOperations, proposal: dict[str, Any], unavailable_starting_materials: list[str] | None,
    include_conditions: bool, include_forward: bool, evidence_refs: list[str] | None,
) -> dict[str, Any]:
    """Retain the submitted hypothesis and full authoritative assessment side by side."""
    refs = _evidence(operations, evidence_refs)
    route = _proposal(proposal)
    unavailable = _strings(unavailable_starting_materials, "unavailable_starting_materials")
    if len(unavailable) > 100:
        raise ValueError("At most 100 unavailable starting materials may be declared")
    assessment = assess_external_route_proposal(
        route, **_assessment_arguments(operations, include_conditions, include_forward),
    )
    return {
        "schema_version": "route_investigation.v1", "origin": "agent_proposal",
        "review_status": "unreviewed", "proposal": route.to_dict(), "assessment": assessment.to_dict(),
        "assessment_options": {"include_conditions": include_conditions, "include_forward": include_forward},
        "material_constraints": assess_declared_route_materials(assessment, tuple(unavailable)).to_dict(),
        "proposed_recipe_assessments": {
            step.external_step_id: _recipe(operations, step.proposal) for step in route.steps
        },
        "evidence_refs": refs, "experimental_feasibility": "not_established",
        "limitations": [
            "An assessed proposal remains a hypothesis, not an observed or experimentally validated route.",
            "Unknown, not_run and out_of_scope checks are not evidence of impossibility or success.",
            "Supplied citations are provenance; their scientific support requires source inspection.",
            "Declared unavailable materials are constraints; actual stock and purchasability are not assessed.",
        ],
    }


def _record(operations: ScientificOperations, source_ref: str) -> dict[str, Any]:
    source = operations.store.read_artifact(source_ref)
    if source.get("operation") not in ROUTE_RECORD_OPERATIONS or source.get("execution_status") != "completed":
        raise ValueError("source_ref must identify a completed route proposal assessment or revision")
    if source["result"].get("schema_version") != "route_investigation.v1":
        raise ValueError("Unsupported route investigation schema")
    return source["result"]


def prepare_planned_route(
    operations: ScientificOperations, source_ref: str, route_id: str,
    unavailable_starting_materials: list[str] | None, include_conditions: bool, include_forward: bool,
) -> dict[str, Any]:
    """Strip trusted planner annotations, then independently assess the selected tree."""
    from core_retrosynthesis.route_contract import ReactionRouteTree

    source = operations.store.read_artifact(source_ref)
    if source.get("operation") not in {"plan_routes", "revise_routes"} or source.get("execution_status") != "completed":
        raise ValueError("source_ref must identify a completed planner call")
    result = source["result"]["response"].get("result") or {}
    routes = [route for route in (*result.get("routes", []), *result.get("partial_routes", [])) if route["route_id"] == route_id]
    if len(routes) != 1:
        raise ValueError("Select exactly one known planner route_id")
    tree = ReactionRouteTree.from_dict(routes[0]["route_tree"])
    proposal = external_route_proposal_from_tree(tree)
    record = assess_route(operations, proposal.to_dict(), unavailable_starting_materials,
                          include_conditions, include_forward, [source_ref])
    return {**record, "source_ref": source_ref, "planner_route_id": route_id}


def inspect_step(operations: ScientificOperations, source_ref: str, step_id: str) -> dict[str, Any]:
    """Retrieve one assessed step, neighboring branches, and canonical molecule audits."""
    record = _record(operations, source_ref)
    route = _proposal(record["proposal"])
    neighbors = external_route_step_neighbors(route, step_id)
    step = next(step for step in route.steps if step.external_step_id == step_id)
    assessment = next(item["assessment"] for item in record["assessment"]["step_assessments"] if item["external_step_id"] == step_id)
    molecules = [assessment["canonical_target_smiles"], *(assessment["canonical_precursor_smiles"] or "").split(".")]
    return {
        "schema_version": "route_step_inspection.v1", "source_ref": source_ref,
        **neighbors, "proposal": step.to_dict(), "assessment": assessment,
        "molecule_audits": {smiles: operations.analyze_molecule(smiles) for smiles in dict.fromkeys(molecules) if smiles},
        "proposed_recipe_assessment": record["proposed_recipe_assessments"].get(step_id),
        "route_topology_gates": record["assessment"]["topology_gates"],
        "material_constraints": record["material_constraints"],
    }


def revise_branch(
    operations: ScientificOperations, source_ref: str, remove_step_ids: list[str],
    replacement_steps: list[dict[str, Any]], reason: str, evidence_refs: list[str] | None,
    assumptions: list[str] | None, risks: list[str] | None,
) -> dict[str, Any]:
    """Record a new proposed branch and recheck the whole route without mutating history."""
    if not isinstance(reason, str) or not reason.strip():
        raise ValueError("A revision requires an explicit reason")
    if not isinstance(replacement_steps, list) or len(replacement_steps) > 40:
        raise ValueError("replacement_steps must be a list of at most 40 physical steps")
    source = _record(operations, source_ref)
    assumptions = _strings(assumptions, "assumptions")
    risks = _strings(risks, "risks")
    if not risks:
        raise ValueError("Revision requires explicit unresolved risks")
    refs = _evidence(operations, [source_ref, *_strings(evidence_refs, "evidence_refs")])
    revision = revise_external_route_proposal(
        _proposal(source["proposal"]), remove_step_ids=tuple(_strings(remove_step_ids, "remove_step_ids")),
        replacement_steps=tuple(_step(value, named=True) for value in replacement_steps),
    )
    record = assess_route(
        operations, revision.proposal.to_dict(), source["material_constraints"]["unavailable_starting_materials"],
        **source["assessment_options"], evidence_refs=refs,
    )
    change = revision.to_dict()
    del change["proposal"]
    return {
        **record, "source_ref": source_ref,
        "revision": {**change, "reason": reason, "assumptions": assumptions, "risks": risks,
                     "basis": "proposed", "reassessment_scope": "all_steps_and_route_topology",
                     "improvement_status": "not_automatically_established"},
    }


def compare_routes(operations: ScientificOperations, source_refs: list[str]) -> dict[str, Any]:
    """Produce an evidence table without collapsing independent checks into a route score."""
    refs = _strings(source_refs, "source_refs")
    if not 2 <= len(refs) <= 5 or len(set(refs)) != len(refs):
        raise ValueError("Compare two to five distinct recorded route alternatives")
    records = [_record(operations, ref) for ref in refs]
    targets = {record["assessment"]["canonical_target_smiles"] for record in records}
    if None in targets or len(targets) != 1:
        raise ValueError("Compared routes must have the same valid canonical target")
    options = records[0]["assessment_options"]
    constraints = records[0]["material_constraints"]["unavailable_starting_materials"]
    if any(record["assessment_options"] != options or
           record["material_constraints"]["unavailable_starting_materials"] != constraints for record in records):
        raise ValueError("Compare routes assessed with the same options and declared material constraints")
    return {
        "schema_version": "route_comparison.v1", "target_smiles": next(iter(targets)),
        "source_refs": refs, "ranking": "not_performed", "experimental_feasibility": "not_established",
        "alternatives": [{
            "source_ref": ref, "route_id": record["assessment"]["route_id"],
            "status": record["assessment"]["status"], "step_count": len(record["proposal"]["steps"]),
            "topology_gates": record["assessment"]["topology_gates"],
            "unresolved_step_ids": record["assessment"]["unresolved_step_ids"],
            "leaf_smiles": record["assessment"]["leaf_smiles"], "material_constraints": record["material_constraints"],
            "assessment_options": record["assessment_options"],
            "steps": [{"step_id": item["external_step_id"],
                       **{key: item["assessment"][key] for key in ("status", "strongest_evidence_tier", "gates", "warnings")}}
                      for item in record["assessment"]["step_assessments"]],
            "proposed_recipe_assessments": record["proposed_recipe_assessments"],
            "revision": record.get("revision"), "limitations": record["limitations"],
        } for ref, record in zip(refs, records)],
    }
