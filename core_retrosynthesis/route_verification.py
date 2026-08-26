"""Independent deterministic verification report for one planned route."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Literal

from .chemistry import canonical_smiles, split_reaction_smiles
from .multistep import MultistepRetrosynthesisRoute
from .route_contract import validate_route_tree
from .route_refinement import collect_route_refinement_issues


ROUTE_VERIFICATION_SCHEMA_VERSION = "route_verification.v1"

RouteVerificationGateStatus = Literal["pass", "warn", "fail"]
RouteVerificationStatus = Literal["verified", "verified_with_cautions", "failed"]


@dataclass(frozen=True)
class RouteVerificationGate:
    """One deterministic whole-route acceptance gate."""

    gate: str
    status: RouteVerificationGateStatus
    message: str
    details: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible gate evidence."""

        return asdict(self)


@dataclass(frozen=True)
class RouteVerificationReport:
    """Auditable whole-route graph and evidence verification."""

    route_id: str
    status: RouteVerificationStatus
    gates: tuple[RouteVerificationGate, ...]
    issue_count: int
    strong_issue_count: int
    schema_version: str = ROUTE_VERIFICATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible verification report."""

        return {
            **asdict(self),
            "gates": [item.to_dict() for item in self.gates],
        }


def _step_consistency_issues(
    route: MultistepRetrosynthesisRoute,
) -> tuple[str, ...]:
    issues: list[str] = []
    for index, step in enumerate(route.steps, start=1):
        product = canonical_smiles(step.product_smiles)
        candidate_product = canonical_smiles(step.candidate.target_smiles)
        precursor = canonical_smiles(".".join(step.precursor_smiles))
        candidate_precursor = canonical_smiles(step.candidate.precursor_smiles)
        reaction = split_reaction_smiles(step.candidate.proposed_reaction_smiles)
        if product is None or product != candidate_product:
            issues.append(f"step-{index}:candidate_product_mismatch")
        if precursor is None or precursor != candidate_precursor:
            issues.append(f"step-{index}:candidate_precursor_mismatch")
        if reaction is None:
            issues.append(f"step-{index}:invalid_proposed_reaction")
            continue
        reaction_precursor = canonical_smiles(reaction[0])
        reaction_product = canonical_smiles(reaction[1])
        if reaction_precursor != precursor:
            issues.append(f"step-{index}:reaction_precursor_mismatch")
        if reaction_product != product:
            issues.append(f"step-{index}:reaction_product_mismatch")
    return tuple(issues)


def verify_planned_route(
    route: MultistepRetrosynthesisRoute,
) -> RouteVerificationReport:
    """Verify graph consistency and evidence gates without changing the route."""

    gates: list[RouteVerificationGate] = []
    tree_report = validate_route_tree(route.route_tree)
    gates.append(
        RouteVerificationGate(
            gate="route_tree_integrity",
            status="pass" if tree_report.valid else "fail",
            message=(
                "The canonical route tree is internally consistent."
                if tree_report.valid
                else "The canonical route tree failed structural validation."
            ),
            details=tuple(tree_report.issues),
        )
    )
    target = canonical_smiles(route.target_smiles)
    tree_target = canonical_smiles(route.route_tree.target_smiles)
    target_valid = target is not None and target == tree_target
    gates.append(
        RouteVerificationGate(
            gate="target_identity",
            status="pass" if target_valid else "fail",
            message=(
                "Route and route-tree targets have the same canonical identity."
                if target_valid
                else "Route and route-tree target identities do not agree."
            ),
        )
    )
    consistency_issues = _step_consistency_issues(route)
    gates.append(
        RouteVerificationGate(
            gate="step_graph_consistency",
            status="pass" if not consistency_issues else "fail",
            message=(
                "Every step agrees with its candidate and reaction graph."
                if not consistency_issues
                else "At least one step disagrees with its stored graph evidence."
            ),
            details=consistency_issues,
        )
    )
    invalid_forward = tuple(
        f"step-{index}:{step.candidate.forward_validation_status}"
        for index, step in enumerate(route.steps, start=1)
        if step.candidate.forward_validation_status != "verified_signature"
    )
    gates.append(
        RouteVerificationGate(
            gate="forward_validation",
            status="pass" if not invalid_forward else "fail",
            message=(
                "Every planned step has verified-signature forward evidence."
                if not invalid_forward
                else "At least one planned step lacks verified forward evidence."
            ),
            details=invalid_forward,
        )
    )
    unresolved = tuple(
        leaf.route_node_id for leaf in route.leaves if not leaf.terminal
    )
    terminal = route.solved and not unresolved
    gates.append(
        RouteVerificationGate(
            gate="terminal_completion",
            status="pass" if terminal else "fail",
            message=(
                "Every route leaf satisfies the terminal predicate."
                if terminal
                else "The route is partial or contains unresolved leaves."
            ),
            details=unresolved,
        )
    )
    issues = collect_route_refinement_issues(route)
    strong = tuple(item.issue_id for item in issues if item.severity == "strong")
    advisory = tuple(item.issue_id for item in issues if item.severity == "advisory")
    issue_status: RouteVerificationGateStatus = (
        "fail" if strong else "warn" if advisory else "pass"
    )
    gates.append(
        RouteVerificationGate(
            gate="chemistry_issue_screen",
            status=issue_status,
            message=(
                "Strong deterministic route issues remain."
                if strong
                else "Only advisory deterministic route issues remain."
                if advisory
                else "No deterministic route issues were detected."
            ),
            details=(*strong, *advisory),
        )
    )
    condition_steps = tuple(
        step for step in route.steps if step.condition_evidence is not None
    )
    condition_gaps = tuple(
        f"step-{index}"
        for index, step in enumerate(route.steps, start=1)
        if step.condition_evidence is not None
        and step.condition_evidence.status == "insufficient_evidence"
    )
    if not condition_steps:
        condition_status: RouteVerificationGateStatus = "warn"
        condition_message = "Conditions were not assessed for the planned steps."
    elif condition_gaps:
        condition_status = "warn"
        condition_message = "Some planned steps lack sufficient condition evidence."
    else:
        condition_status = "pass"
        condition_message = "Every condition-assessed step has supporting evidence."
    gates.append(
        RouteVerificationGate(
            gate="condition_evidence",
            status=condition_status,
            message=condition_message,
            details=condition_gaps,
        )
    )
    if any(item.status == "fail" for item in gates):
        status: RouteVerificationStatus = "failed"
    elif any(item.status == "warn" for item in gates):
        status = "verified_with_cautions"
    else:
        status = "verified"
    return RouteVerificationReport(
        route_id=route.route_id,
        status=status,
        gates=tuple(gates),
        issue_count=len(issues),
        strong_issue_count=len(strong),
    )


__all__ = [
    "ROUTE_VERIFICATION_SCHEMA_VERSION",
    "RouteVerificationGate",
    "RouteVerificationGateStatus",
    "RouteVerificationReport",
    "RouteVerificationStatus",
    "verify_planned_route",
]
