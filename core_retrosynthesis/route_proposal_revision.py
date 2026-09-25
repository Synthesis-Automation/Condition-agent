"""Explicit edits and declared-material checks for review-only route proposals.

These operations do not admit a proposal. Every edited route must pass through
``assess_external_route_proposal`` again, including unchanged downstream steps.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Any

from .chemistry import canonical_smiles
from .external_route_admission import (
    ExternalRouteAssessment, ExternalRouteProposal, ExternalRouteStepProposal,
)


@dataclass(frozen=True)
class ExternalRouteRevision:
    """An explicit replacement retaining the target and every untouched step."""

    proposal: ExternalRouteProposal
    removed_step_ids: tuple[str, ...]
    added_step_ids: tuple[str, ...]
    replaced_step_ids: tuple[str, ...]
    preserved_step_ids: tuple[str, ...]
    schema_version: str = "route_proposal_revision.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the revision without mutable references to its source."""
        return {**asdict(self), "proposal": self.proposal.to_dict()}


def revise_external_route_proposal(
    source: ExternalRouteProposal, *, remove_step_ids: tuple[str, ...],
    replacement_steps: tuple[ExternalRouteStepProposal, ...],
) -> ExternalRouteRevision:
    """Replace explicitly named steps; do not silently prune disconnected branches.

An empty removal set can extend a previously terminal branch. Broken topology
is retained for the canonical route assessor to report, not silently repaired.
"""
    original_ids = tuple(step.external_step_id for step in source.steps)
    replacement_ids = tuple(step.external_step_id for step in replacement_steps)
    for values in (original_ids, remove_step_ids, replacement_ids):
        if len(values) != len(set(values)):
            raise ValueError("Step IDs must be unique within source, removals and replacements")
    if set(remove_step_ids) - set(original_ids):
        raise ValueError("Cannot remove an unknown source step")
    retained = tuple(step for step in source.steps if step.external_step_id not in remove_step_ids)
    preserved_ids = tuple(step.external_step_id for step in retained)
    if set(replacement_ids) & set(preserved_ids):
        raise ValueError("Replacing an existing step requires its explicit removal")
    revised = replace(source, steps=(*retained, *replacement_steps))
    if {step.external_step_id: step.to_dict() for step in source.steps} == {
        step.external_step_id: step.to_dict() for step in revised.steps
    }:
        raise ValueError("Revision must change at least one step")
    return ExternalRouteRevision(
        proposal=revised,
        removed_step_ids=tuple(sorted(set(remove_step_ids) - set(replacement_ids))),
        added_step_ids=tuple(sorted(set(replacement_ids) - set(original_ids))),
        replaced_step_ids=tuple(sorted(set(replacement_ids) & set(remove_step_ids))),
        preserved_step_ids=preserved_ids,
    )


@dataclass(frozen=True)
class DeclaredMaterialAssessment:
    """A check against user-declared unavailable leaves, never stock verification."""

    status: str
    unavailable_starting_materials: tuple[str, ...]
    blocked_leaf_smiles: tuple[str, ...]
    stock_availability: str = "not_assessed"
    schema_version: str = "declared_route_materials.v1"

    def to_dict(self) -> dict[str, Any]:
        """Return the declared constraint scope and matched leaves."""
        return asdict(self)


def assess_declared_route_materials(
    assessment: ExternalRouteAssessment, unavailable_smiles: tuple[str, ...],
) -> DeclaredMaterialAssessment:
    """Check molecular identity of leaves only when route topology is established.

Making an unavailable intermediate in an added step can remove it from the
leaf set. This does not demonstrate availability of the new starting materials.
"""
    normalized = []
    for value in unavailable_smiles:
        molecule = canonical_smiles(value)
        if not molecule or "." in molecule:
            raise ValueError("Each unavailable starting material must be one valid molecular graph")
        normalized.append(molecule)
    declared = tuple(sorted(set(normalized)))
    topology_valid = all(
        gate.status == "pass" for gate in assessment.topology_gates if gate.gate_id != "step_admission"
    )
    parsed = all(
        step.assessment.canonical_target_smiles and step.assessment.canonical_precursor_smiles
        for step in assessment.step_assessments
    )
    if not topology_valid or not parsed:
        return DeclaredMaterialAssessment("unknown", declared, ())
    blocked = tuple(sorted(set(assessment.leaf_smiles) & set(declared)))
    status = "violated" if blocked else "satisfied_for_declared_constraints" if declared else "not_requested"
    return DeclaredMaterialAssessment(status, declared, blocked)


def external_route_step_neighbors(route: ExternalRouteProposal, step_id: str) -> dict[str, Any]:
    """Expose graph-matched producers and consumers; ambiguous matches remain plural."""
    matches = [step for step in route.steps if step.external_step_id == step_id]
    if len(matches) != 1:
        raise ValueError("Select exactly one known external_step_id")
    products = {step.external_step_id: canonical_smiles(step.proposal.target_smiles) for step in route.steps}
    precursors = {
        step.external_step_id: tuple((canonical_smiles(step.proposal.precursor_smiles) or "").split("."))
        for step in route.steps
    }
    return {
        "step_id": step_id,
        "upstream_step_ids": sorted(key for key, product in products.items()
                                    if product and product in precursors[step_id]),
        "downstream_step_ids": sorted(key for key, values in precursors.items()
                                      if products[step_id] and products[step_id] in values),
        "topology_scope": "graph_identity_matches_not_topology_validation",
    }
