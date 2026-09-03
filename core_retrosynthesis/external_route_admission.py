"""Topology-safe review assembly for externally proposed complete routes."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

from forward_synthesis import ForwardOperatorLibrary

from .chemistry import canonical_smiles, digest
from .external_proposal_assessment import (
    ConditionEvidenceEvaluator,
    ExternalProposalAssessmentLimits,
    ExternalProposalGate,
    ExternalProposalSource,
    ExternalRetrosynthesisAssessment,
    ExternalRetrosynthesisProposal,
    assess_external_retrosynthesis_proposal,
    load_external_proposal_admission_policy,
)
from .generic_models import GenericTemplateLibrary
from .route_contract import (
    MoleculeOccurrenceNode,
    PlannedRouteAction,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    assert_valid_route_tree,
    iter_molecule_occurrences,
)


EXTERNAL_ROUTE_ADMISSION_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ExternalRouteStepProposal:
    """One named physical step in an external route proposal."""

    external_step_id: str
    proposal: ExternalRetrosynthesisProposal

    def __post_init__(self) -> None:
        if not self.external_step_id:
            raise ValueError("external route step ID is required")

    def to_dict(self) -> dict[str, Any]:
        return {
            "external_step_id": self.external_step_id,
            **self.proposal.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ExternalRouteStepProposal":
        proposal_value = {
            key: item
            for key, item in value.items()
            if key != "external_step_id"
        }
        return cls(
            external_step_id=str(value.get("external_step_id") or ""),
            proposal=ExternalRetrosynthesisProposal.from_dict(proposal_value),
        )


@dataclass(frozen=True)
class ExternalRouteProposal:
    """A complete external route graph expressed as structured physical steps."""

    target_smiles: str
    steps: tuple[ExternalRouteStepProposal, ...]
    sources: tuple[ExternalProposalSource, ...] = ()
    external_route_id: str | None = None
    schema_version: str = "1.0"

    def to_dict(self) -> dict[str, Any]:
        return {
            "target_smiles": self.target_smiles,
            "steps": [item.to_dict() for item in self.steps],
            "sources": [item.to_dict() for item in self.sources],
            "external_route_id": self.external_route_id,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ExternalRouteProposal":
        if str(value.get("schema_version") or "1.0") != "1.0":
            raise ValueError("unsupported external route proposal schema")
        steps = value.get("steps")
        if not isinstance(steps, list):
            raise ValueError("external route steps must be a list")
        return cls(
            target_smiles=str(value.get("target_smiles") or ""),
            steps=tuple(ExternalRouteStepProposal.from_dict(item) for item in steps),
            sources=tuple(
                ExternalProposalSource.from_dict(item)
                for item in value.get("sources") or ()
            ),
            external_route_id=(
                str(value["external_route_id"])
                if value.get("external_route_id") is not None
                else None
            ),
        )


@dataclass(frozen=True)
class ExternalRouteStepAssessment:
    """External step identity paired with its independent assessment."""

    external_step_id: str
    assessment: ExternalRetrosynthesisAssessment

    def to_dict(self) -> dict[str, Any]:
        return {
            "external_step_id": self.external_step_id,
            "assessment": self.assessment.to_dict(),
        }


@dataclass(frozen=True)
class ExternalRouteAssessment:
    """Review-only route assessment with an optional canonical route tree."""

    assessment_id: str
    route_id: str
    external_route_id: str | None
    canonical_target_smiles: str | None
    status: str
    topology_gates: tuple[ExternalProposalGate, ...]
    step_assessments: tuple[ExternalRouteStepAssessment, ...]
    admitted_route_tree: ReactionRouteTree | None
    admission_eligible: bool
    actionable: bool
    unresolved_step_ids: tuple[str, ...]
    disconnected_step_ids: tuple[str, ...]
    leaf_smiles: tuple[str, ...]
    warnings: tuple[str, ...]
    ranking_influence: str = "none_review_only"
    schema_version: str = EXTERNAL_ROUTE_ADMISSION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            **asdict(self),
            "topology_gates": [item.to_dict() for item in self.topology_gates],
            "step_assessments": [item.to_dict() for item in self.step_assessments],
            "admitted_route_tree": (
                self.admitted_route_tree.to_dict()
                if self.admitted_route_tree is not None
                else None
            ),
        }


def external_route_proposal_from_tree(
    tree: ReactionRouteTree,
    *,
    sources: tuple[ExternalProposalSource, ...] = (),
) -> ExternalRouteProposal:
    """Strip trusted internal identities from a route tree for parity review.

    A retained validated mapped reaction is copied only as correspondence input;
    operator and template annotations are deliberately not copied.
    """

    steps = []
    for molecule in iter_molecule_occurrences(tree):
        reaction = molecule.reaction
        if reaction is None:
            continue
        mapped_reaction = None
        if (
            reaction.evidence.reactants_smiles
            and reaction.evidence.product_smiles_mapped
        ):
            mapped_reaction = (
                f"{reaction.evidence.reactants_smiles}>>"
                f"{reaction.evidence.product_smiles_mapped}"
            )
        steps.append(
            ExternalRouteStepProposal(
                external_step_id=reaction.step_id,
                proposal=ExternalRetrosynthesisProposal(
                    target_smiles=molecule.smiles,
                    precursor_smiles=".".join(
                        child.smiles for child in reaction.children
                    ),
                    mapped_reaction_smiles=mapped_reaction,
                    sources=sources,
                    external_proposal_id=reaction.step_id,
                    claimed_transformation=None,
                ),
            )
        )
    return ExternalRouteProposal(
        target_smiles=tree.target_smiles,
        steps=tuple(steps),
        sources=sources,
        external_route_id=tree.tree_id,
    )


def assess_external_route_proposal(
    route: ExternalRouteProposal,
    operator_library: GenericTemplateLibrary,
    *,
    forward_library: ForwardOperatorLibrary | None = None,
    condition_evaluator: ConditionEvidenceEvaluator | None = None,
    limits: ExternalProposalAssessmentLimits | None = None,
) -> ExternalRouteAssessment:
    """Assess each step independently, then assemble only a consistent route."""

    policy = load_external_proposal_admission_policy()
    bounded = limits or policy.limits
    target = canonical_smiles(route.target_smiles)
    step_ids = tuple(item.external_step_id for item in route.steps)
    if len(route.steps) > bounded.maximum_route_steps:
        raise ValueError(
            f"external route exceeds {bounded.maximum_route_steps} physical steps"
        )
    duplicate_ids = tuple(
        sorted({item for item in step_ids if step_ids.count(item) > 1})
    )
    step_assessments = tuple(
        ExternalRouteStepAssessment(
            external_step_id=item.external_step_id,
            assessment=assess_external_retrosynthesis_proposal(
                item.proposal,
                operator_library,
                forward_library=forward_library,
                condition_evaluator=condition_evaluator,
                limits=bounded,
            ),
        )
        for item in route.steps
    )
    warnings: set[str] = {"EXTERNAL_ROUTE_REVIEW_ONLY"}
    topology_gates: list[ExternalProposalGate] = []
    if target is None or "." in target:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="route_target",
                status="fail",
                summary="Route target must be one valid molecular graph.",
                warnings=("INVALID_EXTERNAL_ROUTE_TARGET",),
            )
        )
    else:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="route_target",
                status="pass",
                summary="Route target parsed as one molecular graph.",
            )
        )
    if not route.steps:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="step_identity",
                status="fail",
                summary="A complete external route requires at least one physical step.",
            )
        )
    elif duplicate_ids:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="step_identity",
                status="fail",
                summary="External step IDs must be unique.",
                evidence_ids=duplicate_ids,
            )
        )
    else:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="step_identity",
                status="pass",
                summary="Every external physical step has a unique occurrence ID.",
            )
        )

    product_to_step: dict[str, ExternalRouteStepAssessment] = {}
    duplicate_products: set[str] = set()
    for step in step_assessments:
        product = step.assessment.canonical_target_smiles
        if product is None:
            continue
        if product in product_to_step:
            duplicate_products.add(product)
        else:
            product_to_step[product] = step
    if duplicate_products:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="product_ownership",
                status="fail",
                summary="More than one external step claims the same product graph.",
                evidence_ids=tuple(sorted(duplicate_products)),
            )
        )
    else:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="product_ownership",
                status="pass",
                summary="Each proposed intermediate product is owned by at most one step.",
            )
        )

    visited: set[str] = set()
    visiting: set[str] = set()
    cycle_steps: set[str] = set()
    leaves: set[str] = set()

    def walk(product: str) -> None:
        step = product_to_step.get(product)
        if step is None:
            leaves.add(product)
            return
        step_id = step.external_step_id
        if step_id in visiting:
            cycle_steps.add(step_id)
            return
        if step_id in visited:
            return
        visiting.add(step_id)
        precursors = step.assessment.canonical_precursor_smiles
        if precursors is not None:
            for precursor in precursors.split("."):
                walk(precursor)
        visiting.remove(step_id)
        visited.add(step_id)

    if target is not None and target in product_to_step and not duplicate_products:
        walk(target)
    disconnected = tuple(sorted(set(step_ids).difference(visited)))
    if target is None or target not in product_to_step:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="route_connectivity",
                status="fail",
                summary="No proposed step produces the declared route target.",
            )
        )
    elif cycle_steps:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="route_connectivity",
                status="fail",
                summary="The proposed route contains a molecular dependency cycle.",
                evidence_ids=tuple(sorted(cycle_steps)),
            )
        )
    elif disconnected:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="route_connectivity",
                status="fail",
                summary="Some proposed steps are disconnected from the target route tree.",
                evidence_ids=disconnected,
            )
        )
    else:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="route_connectivity",
                status="pass",
                summary="Every proposed step is connected to the target in an acyclic route tree.",
            )
        )

    unresolved = tuple(
        step.external_step_id
        for step in step_assessments
        if not step.assessment.admission_eligible
    )
    if unresolved:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="step_admission",
                status="unresolved",
                summary="One or more steps lack the minimum exact-operator and precedent evidence.",
                evidence_ids=unresolved,
            )
        )
    else:
        topology_gates.append(
            ExternalProposalGate(
                gate_id="step_admission",
                status="pass",
                summary="Every proposed physical step meets the review admission threshold.",
            )
        )

    topology_valid = all(item.status != "fail" for item in topology_gates)
    eligible = topology_valid and not unresolved
    route_id = digest(
        "EXTROUTE1",
        target or route.target_smiles,
        *sorted(item.assessment.proposal_id for item in step_assessments),
    )
    tree = (
        _assemble_route_tree(
            route_id,
            target,
            product_to_step,
            sources=route.sources,
        )
        if eligible and target is not None
        else None
    )
    if not topology_valid:
        status = "invalid"
    elif unresolved:
        status = "partially_supported"
    else:
        status = "admitted_review_only"
    assessment_id = digest(
        "EXTROUTEASS1",
        route_id,
        *(item.assessment.assessment_id for item in step_assessments),
        EXTERNAL_ROUTE_ADMISSION_SCHEMA_VERSION,
        policy.definition_id,
    )
    return ExternalRouteAssessment(
        assessment_id=assessment_id,
        route_id=route_id,
        external_route_id=route.external_route_id,
        canonical_target_smiles=target,
        status=status,
        topology_gates=tuple(topology_gates),
        step_assessments=step_assessments,
        admitted_route_tree=tree,
        admission_eligible=eligible,
        actionable=False,
        unresolved_step_ids=unresolved,
        disconnected_step_ids=disconnected,
        leaf_smiles=tuple(sorted(leaves)),
        warnings=tuple(sorted(warnings)),
    )


def _assemble_route_tree(
    route_id: str,
    target: str,
    product_to_step: Mapping[str, ExternalRouteStepAssessment],
    *,
    sources: tuple[ExternalProposalSource, ...],
) -> ReactionRouteTree:
    reaction_count = 0
    maximum_depth = 0
    fingerprint_tokens: set[str] = set()

    def molecule_node(smiles: str, depth: int, path: str) -> MoleculeOccurrenceNode:
        nonlocal reaction_count, maximum_depth
        maximum_depth = max(maximum_depth, depth)
        occurrence_id = digest("EXTMOL1", route_id, path, smiles)
        step = product_to_step.get(smiles)
        if step is None:
            return MoleculeOccurrenceNode(
                occurrence_id=occurrence_id,
                smiles=smiles,
                depth=depth,
                terminal=False,
                terminal_evidence="external_proposal_leaf_stock_not_assessed",
                unresolved_reason="leaf_stock_status_not_assessed",
            )
        assessment = step.assessment
        exact = next(
            item
            for item in assessment.operator_matches
            if item.match_id == assessment.selected_operator_match_id
        )
        assert assessment.canonical_precursor_smiles is not None
        assert assessment.normalized_mapped_reaction_smiles is not None
        assert assessment.reaction_signature is not None
        children = tuple(
            molecule_node(
                precursor,
                depth + 1,
                f"{path}/{step.external_step_id}/{index}",
            )
            for index, precursor in enumerate(
                assessment.canonical_precursor_smiles.split(".")
            )
        )
        reaction_count += 1
        fingerprint_tokens.add(exact.operator_id)
        fingerprint_tokens.add(assessment.proposal_id)
        reaction = RouteReactionNode(
            reaction_node_id=digest(
                "EXTRXN1", route_id, path, step.external_step_id
            ),
            step_id=step.external_step_id,
            depth=depth + 1,
            reaction_smiles=assessment.normalized_mapped_reaction_smiles,
            evidence=RouteStepEvidence(
                evidence_kind="predicted",
                source_dataset_id="external_proposal",
                source_route_id=route_id,
                source_reaction_id=step.external_step_id,
                connectivity_method="external_structures_validated",
                reactants_smiles=assessment.canonical_precursor_smiles,
                product_smiles_mapped=assessment.normalized_mapped_reaction_smiles.split(">>", 1)[1],
                abstracted_reaction_smiles=assessment.normalized_mapped_reaction_smiles,
                abstraction_status="none",
                reaction_signature_id=assessment.reaction_signature.signature_id,
                reaction_signature_schema_version=assessment.reaction_signature.schema_version,
                warnings=("EXTERNAL_PROPOSAL_REVIEW_ONLY",),
            ),
            children=children,
            planned_action=PlannedRouteAction(
                operator_id=exact.operator_id,
                disconnection_site_key=(
                    assessment.reaction_identity.disconnection_site_key
                    if assessment.reaction_identity is not None
                    else ""
                ),
                template_id=assessment.selected_template_id,
            ),
        )
        return MoleculeOccurrenceNode(
            occurrence_id=occurrence_id,
            smiles=smiles,
            depth=depth,
            terminal=False,
            terminal_evidence="none",
            unresolved_reason=None,
            reaction=reaction,
        )

    root = molecule_node(target, 0, "root")
    source_ids = tuple(sorted(item.source_id for item in sources))
    tree = ReactionRouteTree(
        tree_id=digest("EXTTREE1", route_id),
        route_kind="planned",
        target_smiles=target,
        root=root,
        reaction_count=reaction_count,
        maximum_depth=maximum_depth,
        fingerprint_tokens=tuple(sorted(fingerprint_tokens)),
        source_dataset_id="external_proposal",
        source_route_id=route_id,
        source_record_schema_version=EXTERNAL_ROUTE_ADMISSION_SCHEMA_VERSION,
        connectivity_method="external_structures_validated",
        warnings=(
            "EXTERNAL_ROUTE_REVIEW_ONLY",
            "LEAF_STOCK_STATUS_NOT_ASSESSED",
            *(f"SOURCE:{source_id}" for source_id in source_ids),
        ),
    )
    assert_valid_route_tree(tree)
    return tree


__all__ = [
    "EXTERNAL_ROUTE_ADMISSION_SCHEMA_VERSION",
    "ExternalRouteAssessment",
    "ExternalRouteProposal",
    "ExternalRouteStepAssessment",
    "ExternalRouteStepProposal",
    "assess_external_route_proposal",
    "external_route_proposal_from_tree",
]
