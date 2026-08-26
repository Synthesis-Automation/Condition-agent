"""Typed, deterministic contracts for bounded route refinement.

The language-model layer may select an objective and a closed refinement
method, but it never supplies molecular structures or graph edits.  Concrete
candidate exclusions are derived here from an already validated route step.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from typing import Any, Literal, Sequence, TYPE_CHECKING

from reactive_taxonomy import assess_reaction_compatibility

from .generic_models import GenericDisconnectionCandidate

if TYPE_CHECKING:
    from .multistep import (
        MultistepRetrosynthesisResult,
        MultistepRetrosynthesisRoute,
    )


ROUTE_REFINEMENT_SCHEMA_VERSION = "route_refinement.v3"

RouteIssueKind = Literal[
    "precursor_compatibility",
    "reaction_compatibility",
    "selectivity",
    "condition_gap",
    "unresolved_leaf",
    "tactical_step",
]
RouteIssueSeverity = Literal["advisory", "strong"]
RouteRefinementObjective = Literal[
    "resolve_compatibility_conflict",
    "resolve_selectivity_warning",
    "resolve_condition_gap",
    "resolve_unresolved_leaf",
    "reduce_tactical_churn",
    "find_alternative_route",
]
RouteRefinementMethod = Literal[
    "alternate_disconnection",
    "alternate_realization",
]
RouteRefinementStatus = Literal[
    "improved_alternative_found",
    "alternatives_found_no_verified_improvement",
    "no_alternative_found",
]
RouteRepairKind = Literal[
    "alternate_realization",
    "alternate_disconnection",
    "temporary_masking_sequence",
]
RouteRepairProposalStatus = Literal["actionable", "unavailable"]

ROUTE_REFINEMENT_OBJECTIVES = frozenset(
    {
        "resolve_compatibility_conflict",
        "resolve_selectivity_warning",
        "resolve_condition_gap",
        "resolve_unresolved_leaf",
        "reduce_tactical_churn",
        "find_alternative_route",
    }
)
ROUTE_REFINEMENT_METHODS = frozenset(
    {"alternate_disconnection", "alternate_realization"}
)

_OBJECTIVE_ISSUE_KINDS: dict[str, frozenset[str]] = {
    "resolve_compatibility_conflict": frozenset(
        {"precursor_compatibility", "reaction_compatibility"}
    ),
    "resolve_selectivity_warning": frozenset({"selectivity"}),
    "resolve_condition_gap": frozenset({"condition_gap"}),
    "resolve_unresolved_leaf": frozenset({"unresolved_leaf"}),
    "reduce_tactical_churn": frozenset({"tactical_step"}),
    "find_alternative_route": frozenset(),
}


def _stable_id(prefix: str, value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    return f"{prefix}:{hashlib.sha256(encoded).hexdigest()[:20]}"


@dataclass(frozen=True)
class RouteRefinementIssue:
    """One graph- or evidence-derived route issue exposed to the controller."""

    issue_id: str
    kind: RouteIssueKind
    severity: RouteIssueSeverity
    subject_type: Literal["route", "step", "leaf"]
    subject_id: str
    message: str
    step_index: int | None = None
    details: tuple[tuple[str, Any], ...] = ()
    schema_version: str = ROUTE_REFINEMENT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.issue_id or not self.subject_id or not self.message.strip():
            raise ValueError("route-refinement issues require IDs and a message")
        if self.kind not in {
            "precursor_compatibility",
            "reaction_compatibility",
            "selectivity",
            "condition_gap",
            "unresolved_leaf",
            "tactical_step",
        }:
            raise ValueError(f"unsupported route issue kind: {self.kind}")
        if self.severity not in {"advisory", "strong"}:
            raise ValueError(f"unsupported route issue severity: {self.severity}")
        if self.subject_type not in {"route", "step", "leaf"}:
            raise ValueError(f"unsupported route issue subject: {self.subject_type}")
        if self.step_index is not None and self.step_index < 1:
            raise ValueError("route issue step index must be positive")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible issue without losing typed evidence."""

        return {
            **asdict(self),
            "details": dict(self.details),
        }


@dataclass(frozen=True)
class RouteRepairProposal:
    """One deterministic action option derived from a typed route issue."""

    proposal_id: str
    issue_id: str
    source_route_id: str
    source_step_id: str
    repair_kind: RouteRepairKind
    status: RouteRepairProposalStatus
    objective: RouteRefinementObjective
    refinement_method: RouteRefinementMethod | None
    maximum_added_steps: int
    reason: str
    schema_version: str = ROUTE_REFINEMENT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not all(
            (
                self.proposal_id,
                self.issue_id,
                self.source_route_id,
                self.source_step_id,
                self.reason.strip(),
            )
        ):
            raise ValueError("route repair proposals require stable provenance")
        if self.status == "actionable" and self.refinement_method is None:
            raise ValueError("actionable repair proposals require a method")
        if self.status == "unavailable" and self.refinement_method is not None:
            raise ValueError("unavailable repair proposals cannot expose a method")
        if self.maximum_added_steps < 0:
            raise ValueError("repair proposal added steps must be nonnegative")

    def to_dict(self) -> dict[str, Any]:
        """Return a structure-free repair proposal."""

        return asdict(self)


@dataclass(frozen=True)
class RouteCandidateExclusion:
    """One internally derived candidate exclusion for a product expansion."""

    product_smiles: str
    strategy_id: str
    realization_id: str | None = None

    def __post_init__(self) -> None:
        if not self.product_smiles or not self.strategy_id:
            raise ValueError("candidate exclusions require product and strategy IDs")

    def matches(
        self,
        product_smiles: str,
        candidate: GenericDisconnectionCandidate,
    ) -> bool:
        """Return whether the validated candidate is excluded by this rule."""

        if product_smiles != self.product_smiles:
            return False
        if candidate.strategy_id != self.strategy_id:
            return False
        return self.realization_id is None or (
            candidate.realization_id == self.realization_id
        )


@dataclass(frozen=True)
class RouteRefinementIntent:
    """A bounded model-selected objective over existing deterministic IDs."""

    source_route_id: str
    source_step_id: str
    objective: RouteRefinementObjective
    method: RouteRefinementMethod
    issue_ids: tuple[str, ...]
    maximum_added_steps: int = 0
    schema_version: str = ROUTE_REFINEMENT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.source_route_id or not self.source_step_id:
            raise ValueError("route refinement requires source route and step IDs")
        if self.objective not in ROUTE_REFINEMENT_OBJECTIVES:
            raise ValueError(f"unsupported route refinement objective: {self.objective}")
        if self.method not in ROUTE_REFINEMENT_METHODS:
            raise ValueError(f"unsupported route refinement method: {self.method}")
        if not self.issue_ids:
            raise ValueError("route refinement must cite at least one typed issue")
        if len(self.issue_ids) != len(set(self.issue_ids)):
            raise ValueError("route refinement issue IDs must be unique")
        if self.maximum_added_steps not in {0, 1}:
            raise ValueError("maximum_added_steps must be zero or one")

    @property
    def intent_id(self) -> str:
        """Return a deterministic identity for replay and lineage."""

        return _stable_id(
            "RINT1",
            {
                "source_route_id": self.source_route_id,
                "source_step_id": self.source_step_id,
                "objective": self.objective,
                "method": self.method,
                "issue_ids": self.issue_ids,
                "maximum_added_steps": self.maximum_added_steps,
                "schema_version": self.schema_version,
            },
        )

    def to_dict(self) -> dict[str, Any]:
        """Return the structure-free controller intent."""

        return {"intent_id": self.intent_id, **asdict(self)}


def enumerate_route_repair_proposals(
    route: "MultistepRetrosynthesisRoute",
    issue: RouteRefinementIssue,
) -> tuple[RouteRepairProposal, ...]:
    """Enumerate bounded repair actions without generating molecular structures."""

    objective_by_kind: dict[str, RouteRefinementObjective] = {
        "precursor_compatibility": "resolve_compatibility_conflict",
        "reaction_compatibility": "resolve_compatibility_conflict",
        "selectivity": "resolve_selectivity_warning",
        "condition_gap": "resolve_condition_gap",
        "unresolved_leaf": "resolve_unresolved_leaf",
        "tactical_step": "reduce_tactical_churn",
    }
    objective = objective_by_kind[issue.kind]
    step = next(
        (value for value in route.steps if value.step_id == issue.subject_id),
        None,
    )
    if step is None:
        return ()

    def proposal(
        repair_kind: RouteRepairKind,
        status: RouteRepairProposalStatus,
        method: RouteRefinementMethod | None,
        maximum_added_steps: int,
        reason: str,
    ) -> RouteRepairProposal:
        identity = {
            "route_id": route.route_id,
            "step_id": step.step_id,
            "issue_id": issue.issue_id,
            "repair_kind": repair_kind,
            "status": status,
            "objective": objective,
            "method": method,
            "maximum_added_steps": maximum_added_steps,
            "schema_version": ROUTE_REFINEMENT_SCHEMA_VERSION,
        }
        return RouteRepairProposal(
            proposal_id=_stable_id("RPROP1", identity),
            issue_id=issue.issue_id,
            source_route_id=route.route_id,
            source_step_id=step.step_id,
            repair_kind=repair_kind,
            status=status,
            objective=objective,
            refinement_method=method,
            maximum_added_steps=maximum_added_steps,
            reason=reason,
        )

    proposals: list[RouteRepairProposal] = []
    if getattr(step.candidate, "realization_id", ""):
        proposals.append(
            proposal(
                "alternate_realization",
                "actionable",
                "alternate_realization",
                1,
                "Search validated concrete realizations while excluding only "
                "the source realization.",
            )
        )
    if getattr(step.candidate, "strategy_id", ""):
        proposals.append(
            proposal(
                "alternate_disconnection",
                "actionable",
                "alternate_disconnection",
                1,
                "Search validated alternatives while excluding the source "
                "disconnection strategy.",
            )
        )
    if issue.kind == "reaction_compatibility":
        rule_ids = tuple(dict(issue.details).get("rule_ids") or ())
        if "RCI001_HALO_CARBON_TRANSFER_PROTIC_CONFLICT" in rule_ids:
            proposals.append(
                proposal(
                    "temporary_masking_sequence",
                    "unavailable",
                    None,
                    2,
                    "No registered protection/deprotection sequence operator is "
                    "available; the system will not invent protected structures.",
                )
            )
    return tuple(proposals)


@dataclass(frozen=True)
class RouteRefinementPlan:
    """Concrete deterministic search policy derived from a route intent."""

    intent: RouteRefinementIntent
    exclusion: RouteCandidateExclusion
    excluded_candidate_scope: Literal["strategy", "realization"]
    schema_version: str = ROUTE_REFINEMENT_SCHEMA_VERSION


@dataclass(frozen=True)
class RouteRefinementCandidateAssessment:
    """Objective comparison for one regenerated route candidate."""

    route_id: str
    solved: bool
    route_cost: float
    reaction_count: int
    relevant_issue_count: int
    source_relevant_issue_count: int
    strong_issue_count: int
    source_strong_issue_count: int
    total_issue_count: int
    source_total_issue_count: int
    verification_status: Literal[
        "verified", "verified_with_cautions", "failed"
    ]
    hard_verification_failures: tuple[str, ...]
    previously_accepted: bool
    improved: bool
    rejection_reasons: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible comparison."""

        return asdict(self)


@dataclass(frozen=True)
class RouteRefinementOutcome:
    """Immutable lineage and deterministic comparison for one re-search."""

    intent: RouteRefinementIntent
    status: RouteRefinementStatus
    source_route_id: str
    candidate_assessments: tuple[RouteRefinementCandidateAssessment, ...]
    accepted_route_id: str | None = None
    retained_route_id: str = ""
    source_route_preserved: bool = True
    schema_version: str = ROUTE_REFINEMENT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return structure-free refinement evidence for the LLM controller."""

        return {
            "schema_version": self.schema_version,
            "intent": self.intent.to_dict(),
            "status": self.status,
            "source_route_id": self.source_route_id,
            "accepted_route_id": self.accepted_route_id,
            "retained_route_id": self.retained_route_id,
            "source_route_preserved": self.source_route_preserved,
            "candidate_assessments": [
                item.to_dict() for item in self.candidate_assessments
            ],
        }


def collect_route_refinement_issues(
    route: "MultistepRetrosynthesisRoute",
) -> tuple[RouteRefinementIssue, ...]:
    """Derive typed route issues without reaction-name or model inference."""

    issues: list[RouteRefinementIssue] = []
    for step_index, step in enumerate(route.steps, start=1):
        candidate = step.candidate
        disposition = candidate.precursor_compatibility_disposition
        if disposition != "pass":
            payload = {
                "route_id": route.route_id,
                "step_id": step.step_id,
                "kind": "precursor_compatibility",
                "disposition": disposition,
                "assessment_ids": tuple(
                    item.assessment_id
                    for item in candidate.precursor_compatibility_assessments
                ),
            }
            issues.append(
                RouteRefinementIssue(
                    issue_id=_stable_id("RISS1", payload),
                    kind="precursor_compatibility",
                    severity=(
                        "strong"
                        if candidate.precursor_compatibility_warning_strength
                        == "strong"
                        else "advisory"
                    ),
                    subject_type="step",
                    subject_id=step.step_id,
                    step_index=step_index,
                    message=(
                        "The deterministic precursor compatibility policy returned "
                        f"{disposition}."
                    ),
                    details=(
                        ("disposition", disposition),
                        (
                            "assessment_ids",
                            payload["assessment_ids"],
                        ),
                        (
                            "policy_definition_id",
                            candidate.precursor_compatibility_policy_definition_id,
                        ),
                    ),
                )
            )
        reaction_smiles = str(
            getattr(candidate, "condition_query_reaction_smiles", "")
            or getattr(candidate, "proposed_reaction_smiles", "")
        )
        stored_reaction_policy = str(
            getattr(
                candidate,
                "reaction_compatibility_policy_definition_id",
                "",
            )
        )
        reaction_assessments = tuple(
            getattr(candidate, "reaction_compatibility_assessments", ())
        )
        if reaction_smiles and not stored_reaction_policy:
            reaction_assessments = assess_reaction_compatibility(reaction_smiles)
        if reaction_assessments:
            assessment_ids = tuple(
                item.assessment_id for item in reaction_assessments
            )
            rule_ids = tuple(
                sorted({item.rule_id for item in reaction_assessments})
            )
            regimes = tuple(
                sorted({item.inferred_regime for item in reaction_assessments})
            )
            payload = {
                "route_id": route.route_id,
                "step_id": step.step_id,
                "kind": "reaction_compatibility",
                "assessment_ids": assessment_ids,
            }
            issues.append(
                RouteRefinementIssue(
                    issue_id=_stable_id("RISS1", payload),
                    kind="reaction_compatibility",
                    severity=(
                        "strong"
                        if any(
                            item.warning_strength == "strong"
                            for item in reaction_assessments
                        )
                        else "advisory"
                    ),
                    subject_type="step",
                    subject_id=step.step_id,
                    step_index=step_index,
                    message=(
                        reaction_assessments[0].message
                        if len(reaction_assessments) == 1
                        else "The deterministic reaction-regime analysis returned "
                        "compatibility conflicts."
                    ),
                    details=(
                        ("assessment_ids", assessment_ids),
                        ("rule_ids", rule_ids),
                        ("inferred_regimes", regimes),
                        (
                            "definition_id",
                            reaction_assessments[0].definition_id,
                        ),
                    ),
                )
            )
        if candidate.selectivity_warnings:
            warning_ids = tuple(
                f"{item.audit_id}:{item.code}"
                for item in candidate.selectivity_warnings
            )
            payload = {
                "route_id": route.route_id,
                "step_id": step.step_id,
                "kind": "selectivity",
                "warning_ids": warning_ids,
            }
            issues.append(
                RouteRefinementIssue(
                    issue_id=_stable_id("RISS1", payload),
                    kind="selectivity",
                    severity="advisory",
                    subject_type="step",
                    subject_id=step.step_id,
                    step_index=step_index,
                    message="The deterministic step analysis returned selectivity warnings.",
                    details=(("warning_ids", warning_ids),),
                )
            )
        condition = step.condition_evidence
        if condition is not None and condition.status == "insufficient_evidence":
            payload = {
                "route_id": route.route_id,
                "step_id": step.step_id,
                "kind": "condition_gap",
                "status": condition.status,
            }
            issues.append(
                RouteRefinementIssue(
                    issue_id=_stable_id("RISS1", payload),
                    kind="condition_gap",
                    severity="advisory",
                    subject_type="step",
                    subject_id=step.step_id,
                    step_index=step_index,
                    message="No sufficient deterministic condition evidence was retrieved.",
                    details=(
                        ("condition_status", condition.status),
                        ("condition_warnings", condition.warnings),
                    ),
                )
            )
        if candidate.strategic_class in {
            "functional_group_interconversion",
            "protection_state_change",
        }:
            payload = {
                "route_id": route.route_id,
                "step_id": step.step_id,
                "kind": "tactical_step",
                "strategic_class": candidate.strategic_class,
            }
            issues.append(
                RouteRefinementIssue(
                    issue_id=_stable_id("RISS1", payload),
                    kind="tactical_step",
                    severity="advisory",
                    subject_type="step",
                    subject_id=step.step_id,
                    step_index=step_index,
                    message=(
                        "The deterministic strategic-complexity analysis classified "
                        "this as a tactical state-change step."
                    ),
                    details=(("strategic_class", candidate.strategic_class),),
                )
            )
    for leaf in route.leaves:
        if leaf.terminal:
            continue
        payload = {
            "route_id": route.route_id,
            "leaf_id": leaf.route_node_id,
            "kind": "unresolved_leaf",
            "reason": leaf.unresolved_reason,
        }
        issues.append(
            RouteRefinementIssue(
                issue_id=_stable_id("RISS1", payload),
                kind="unresolved_leaf",
                severity="strong",
                subject_type="leaf",
                subject_id=leaf.route_node_id,
                message="The route leaf does not satisfy the terminal predicate.",
                details=(("unresolved_reason", leaf.unresolved_reason),),
            )
        )
    return tuple(issues)


def build_route_refinement_plan(
    route: "MultistepRetrosynthesisRoute",
    intent: RouteRefinementIntent,
    issues: Sequence[RouteRefinementIssue],
) -> RouteRefinementPlan:
    """Resolve a structure-free intent to one exact deterministic exclusion."""

    if route.route_id != intent.source_route_id:
        raise ValueError("refinement intent references a different source route")
    issue_by_id = {item.issue_id: item for item in issues}
    unknown = set(intent.issue_ids) - set(issue_by_id)
    if unknown:
        raise ValueError(f"refinement intent cites unknown issue IDs: {sorted(unknown)}")
    expected_kinds = _OBJECTIVE_ISSUE_KINDS[intent.objective]
    cited_kinds = {issue_by_id[value].kind for value in intent.issue_ids}
    if expected_kinds and not expected_kinds.intersection(cited_kinds):
        raise ValueError(
            "refinement objective is not supported by the cited issue kinds"
        )
    step = next(
        (item for item in route.steps if item.step_id == intent.source_step_id),
        None,
    )
    if step is None:
        raise ValueError("refinement intent references an unknown source step")
    for issue_id in intent.issue_ids:
        issue = issue_by_id[issue_id]
        if issue.step_index is not None and issue.subject_id != step.step_id:
            raise ValueError("step-scoped refinement issue references another step")
    candidate = step.candidate
    if not candidate.strategy_id:
        raise ValueError("source candidate has no stable strategy identity")
    realization_id: str | None = None
    scope: Literal["strategy", "realization"] = "strategy"
    if intent.method == "alternate_realization":
        if not candidate.realization_id:
            raise ValueError("source candidate has no stable realization identity")
        realization_id = candidate.realization_id
        scope = "realization"
    return RouteRefinementPlan(
        intent=intent,
        exclusion=RouteCandidateExclusion(
            product_smiles=step.product_smiles,
            strategy_id=candidate.strategy_id,
            realization_id=realization_id,
        ),
        excluded_candidate_scope=scope,
    )


def summarize_route_refinement(
    source_route: "MultistepRetrosynthesisRoute",
    intent: RouteRefinementIntent,
    result: "MultistepRetrosynthesisResult",
    *,
    previously_accepted_route_ids: Sequence[str] = (),
) -> RouteRefinementOutcome:
    """Conservatively compare regenerated routes and select one verified improvement.

    Repair is objective-specific: a candidate must reduce the selected issue kind,
    pass every hard graph gate, preserve solved status, and avoid increasing either
    strong or total deterministic issues.  Previously accepted route identities are
    rejected to keep a multi-attempt controller from cycling.
    """

    # Imported lazily because route verification consumes the issue collector in
    # this module.  The dependency remains entirely inside core_retrosynthesis.
    from .route_verification import verify_planned_route

    relevant_kinds = _OBJECTIVE_ISSUE_KINDS[intent.objective]
    source_issues = collect_route_refinement_issues(source_route)
    source_count = sum(
        not relevant_kinds or item.kind in relevant_kinds for item in source_issues
    )
    source_strong_count = sum(item.severity == "strong" for item in source_issues)
    source_total_count = len(source_issues)
    accepted_history = frozenset(previously_accepted_route_ids)
    candidates = (*result.routes, *result.partial_routes)
    assessments = []
    for route in candidates:
        issues = collect_route_refinement_issues(route)
        issue_count = sum(
            not relevant_kinds or item.kind in relevant_kinds for item in issues
        )
        strong_count = sum(item.severity == "strong" for item in issues)
        report = verify_planned_route(route)
        hard_failures = tuple(
            gate.gate
            for gate in report.gates
            if gate.gate
            in {
                "route_tree_integrity",
                "target_identity",
                "step_graph_consistency",
                "forward_validation",
            }
            and gate.status == "fail"
        )
        previously_accepted = route.route_id in accepted_history
        rejection_reasons: list[str] = []
        if route.route_id == source_route.route_id:
            rejection_reasons.append("source_route_repeated")
        if previously_accepted:
            rejection_reasons.append("previously_accepted_route_repeated")
        if hard_failures:
            rejection_reasons.append("hard_verification_failed")
        if source_route.solved and not route.solved:
            rejection_reasons.append("solved_status_regressed")
        if strong_count > source_strong_count:
            rejection_reasons.append("strong_issue_count_increased")
        if len(issues) > source_total_count:
            rejection_reasons.append("total_issue_count_increased")
        if relevant_kinds and issue_count >= source_count:
            rejection_reasons.append("objective_issue_count_not_reduced")
        improved = not rejection_reasons
        assessments.append(
            RouteRefinementCandidateAssessment(
                route_id=route.route_id,
                solved=route.solved,
                route_cost=route.route_cost,
                reaction_count=route.reaction_count,
                relevant_issue_count=issue_count,
                source_relevant_issue_count=source_count,
                strong_issue_count=strong_count,
                source_strong_issue_count=source_strong_count,
                total_issue_count=len(issues),
                source_total_issue_count=source_total_count,
                verification_status=report.status,
                hard_verification_failures=hard_failures,
                previously_accepted=previously_accepted,
                improved=improved,
                rejection_reasons=tuple(rejection_reasons),
            )
        )
    improving = [item for item in assessments if item.improved]
    accepted = (
        min(
            improving,
            key=lambda item: (
                not item.solved,
                item.strong_issue_count,
                item.relevant_issue_count,
                item.total_issue_count,
                item.reaction_count,
                item.route_cost,
                item.route_id,
            ),
        )
        if improving
        else None
    )
    if accepted is not None:
        status: RouteRefinementStatus = "improved_alternative_found"
    elif assessments:
        status = "alternatives_found_no_verified_improvement"
    else:
        status = "no_alternative_found"
    return RouteRefinementOutcome(
        intent=intent,
        status=status,
        source_route_id=source_route.route_id,
        candidate_assessments=tuple(assessments),
        accepted_route_id=accepted.route_id if accepted is not None else None,
        retained_route_id=(
            accepted.route_id if accepted is not None else source_route.route_id
        ),
    )


__all__ = [
    "ROUTE_REFINEMENT_METHODS",
    "ROUTE_REFINEMENT_OBJECTIVES",
    "ROUTE_REFINEMENT_SCHEMA_VERSION",
    "RouteCandidateExclusion",
    "RouteRepairKind",
    "RouteRepairProposal",
    "RouteRepairProposalStatus",
    "RouteRefinementCandidateAssessment",
    "RouteRefinementIntent",
    "RouteRefinementIssue",
    "RouteRefinementMethod",
    "RouteRefinementObjective",
    "RouteRefinementOutcome",
    "RouteRefinementPlan",
    "build_route_refinement_plan",
    "collect_route_refinement_issues",
    "enumerate_route_repair_proposals",
    "summarize_route_refinement",
]
