"""Advisory forward-competition audits for retrosynthesis candidates and routes."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from typing import Any, Mapping, Tuple

from forward_synthesis import (
    ForwardOperatorLibrary,
    RouteStepForwardAssessment,
    assess_proposed_step,
    build_forward_library,
)

from .generic_models import GenericDisconnectionCandidate, GenericTemplateLibrary
from .route_contract import ReactionRouteTree, iter_molecule_occurrences


ROUTE_FORWARD_ASSESSMENT_SCHEMA_VERSION = "1.1"


@dataclass(frozen=True)
class RouteForwardStepAssessment:
    """One route reaction occurrence and its independent forward audit."""

    step_id: str
    reaction_node_id: str
    product_occurrence_id: str
    assessment: RouteStepForwardAssessment

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RouteForwardAssessmentReport:
    """Advisory forward-competition evidence for every route step."""

    tree_id: str
    step_assessments: Tuple[RouteForwardStepAssessment, ...]
    disposition_counts: Mapping[str, int]
    high_risk_step_ids: Tuple[str, ...]
    route_ranking_impact: str = "none_advisory_only"
    schema_version: str = ROUTE_FORWARD_ASSESSMENT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def assess_retrosynthesis_candidate_forward(
    candidate: GenericDisconnectionCandidate,
    forward_library: ForwardOperatorLibrary,
    *,
    recipe: Mapping[str, Any] | None = None,
    top_k: int = 20,
) -> RouteStepForwardAssessment:
    """Independently challenge one structurally validated disconnection."""

    return assess_proposed_step(
        candidate.precursor_smiles,
        candidate.target_smiles,
        forward_library,
        recipe=recipe,
        operator_hint=candidate.operator_id or None,
        top_k=top_k,
    )


def build_forward_library_from_generic(
    library: GenericTemplateLibrary,
    *,
    require_source_round_trip: bool = True,
) -> ForwardOperatorLibrary:
    """Create a forward-admitted view of the generic operator library."""

    return build_forward_library(
        library,
        require_source_round_trip=require_source_round_trip,
    )


def assess_route_tree_forward(
    tree: ReactionRouteTree,
    forward_library: ForwardOperatorLibrary,
    *,
    recipes_by_step_id: Mapping[str, Mapping[str, Any]] | None = None,
    top_k: int = 20,
) -> RouteForwardAssessmentReport:
    """Run target-blind product competition for every concrete route step."""

    recipes = recipes_by_step_id or {}
    assessments = []
    for product_node in iter_molecule_occurrences(tree):
        reaction = product_node.reaction
        if reaction is None:
            continue
        starting_materials = ".".join(
            sorted(child.smiles for child in reaction.children)
        )
        operator_hint = (
            reaction.planned_action.operator_id
            if reaction.planned_action is not None
            else None
        )
        assessment = assess_proposed_step(
            starting_materials,
            product_node.smiles,
            forward_library,
            recipe=recipes.get(reaction.step_id),
            operator_hint=operator_hint,
            top_k=top_k,
        )
        assessments.append(
            RouteForwardStepAssessment(
                step_id=reaction.step_id,
                reaction_node_id=reaction.reaction_node_id,
                product_occurrence_id=product_node.occurrence_id,
                assessment=assessment,
            )
        )
    counts = Counter(item.assessment.disposition for item in assessments)
    high_risk = tuple(
        item.step_id
        for item in assessments
        if item.assessment.disposition
        in {"competitive", "structurally_inconsistent", "condition_incompatible"}
    )
    return RouteForwardAssessmentReport(
        tree_id=tree.tree_id,
        step_assessments=tuple(assessments),
        disposition_counts=dict(sorted(counts.items())),
        high_risk_step_ids=high_risk,
    )


__all__ = [
    "ROUTE_FORWARD_ASSESSMENT_SCHEMA_VERSION",
    "RouteForwardAssessmentReport",
    "RouteForwardStepAssessment",
    "assess_retrosynthesis_candidate_forward",
    "assess_route_tree_forward",
    "build_forward_library_from_generic",
]
