"""Deterministic issue, repair, and verification contracts for one-step retro.

The module never generates a molecule, reaction, or condition.  It turns
already validated strategies, retained realizations, and canonical condition
recommendations into stable repair choices that an advisory controller may
select by ID.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from types import SimpleNamespace
from typing import Any, Literal, Mapping, Sequence

from .chemistry import canonical_smiles, split_reaction_smiles
from .condition_selectivity_repair import (
    ConditionSelectivityRepairAssessment,
    assess_condition_selectivity_repairs,
)
from .generic_models import GenericDisconnectionCandidate, StrategyProposal
from .route_refinement import (
    RouteRefinementIssue,
    collect_route_refinement_issues,
)


SINGLE_STEP_REFINEMENT_SCHEMA_VERSION = "single_step_refinement.v1"

SingleStepRepairKind = Literal[
    "alternate_realization",
    "alternate_strategy",
    "condition_selectivity",
]
SingleStepRepairStatus = Literal["actionable", "unavailable"]


def _stable_id(prefix: str, value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    return f"{prefix}:{hashlib.sha256(encoded).hexdigest()[:20]}"


def _reaction_identity(candidate: GenericDisconnectionCandidate) -> str:
    return (
        candidate.condition_query_reaction_smiles
        or candidate.proposed_reaction_smiles
    )


def _synthetic_route(
    strategy: StrategyProposal,
    *,
    candidate: GenericDisconnectionCandidate | None = None,
    condition_evidence: Any = None,
    condition_selectivity_assessment: ConditionSelectivityRepairAssessment
    | None = None,
) -> Any:
    selected = candidate or strategy.representative
    step_id = selected.realization_id or strategy.strategy_id
    step = SimpleNamespace(
        step_id=step_id,
        candidate=selected,
        condition_evidence=condition_evidence,
        condition_selectivity_assessment=condition_selectivity_assessment,
    )
    return SimpleNamespace(
        route_id=f"single-step:{strategy.strategy_id}:{step_id}",
        steps=(step,),
        leaves=(),
    )


def collect_single_step_refinement_issues(
    strategy: StrategyProposal,
    *,
    candidate: GenericDisconnectionCandidate | None = None,
    condition_evidence: Any = None,
    condition_selectivity_assessment: ConditionSelectivityRepairAssessment
    | None = None,
) -> tuple[RouteRefinementIssue, ...]:
    """Reuse the canonical step chemistry audit for one retained strategy."""

    return collect_route_refinement_issues(
        _synthetic_route(
            strategy,
            candidate=candidate,
            condition_evidence=condition_evidence,
            condition_selectivity_assessment=condition_selectivity_assessment,
        )
    )


@dataclass(frozen=True)
class SingleStepRepairProposal:
    """One structure-free choice among existing deterministic artifacts."""

    proposal_id: str
    issue_id: str
    source_strategy_id: str
    source_realization_id: str
    repair_kind: SingleStepRepairKind
    status: SingleStepRepairStatus
    objective: str
    reason: str
    target_strategy_id: str | None = None
    target_realization_id: str | None = None
    assessment_id: str | None = None
    details: tuple[tuple[str, Any], ...] = ()
    schema_version: str = SINGLE_STEP_REFINEMENT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not all(
            (
                self.proposal_id,
                self.issue_id,
                self.source_strategy_id,
                self.source_realization_id,
                self.reason.strip(),
            )
        ):
            raise ValueError("single-step repair proposal provenance is incomplete")
        if self.status == "actionable" and not (
            self.target_strategy_id or self.assessment_id
        ):
            raise ValueError("actionable repair proposal requires a deterministic target")

    def to_dict(self) -> dict[str, Any]:
        return {**asdict(self), "details": dict(self.details)}


_OBJECTIVE_BY_ISSUE = {
    "precursor_compatibility": "resolve_compatibility_conflict",
    "reaction_compatibility": "resolve_compatibility_conflict",
    "selectivity": "resolve_selectivity_warning",
    "condition_gap": "resolve_condition_gap",
    "tactical_step": "reduce_tactical_churn",
}


def _issue_profile(issues: Sequence[RouteRefinementIssue], kind: str) -> tuple[int, int, int]:
    relevant = sum(item.kind == kind for item in issues)
    strong = sum(item.severity == "strong" for item in issues)
    return relevant, strong, len(issues)


def _is_improvement(
    source: Sequence[RouteRefinementIssue],
    target: Sequence[RouteRefinementIssue],
    kind: str,
) -> bool:
    source_relevant, source_strong, source_total = _issue_profile(source, kind)
    target_relevant, target_strong, target_total = _issue_profile(target, kind)
    return (
        target_relevant < source_relevant
        and target_strong <= source_strong
        and target_total <= source_total
    )


def enumerate_single_step_repair_proposals(
    source_strategy: StrategyProposal,
    issue: RouteRefinementIssue,
    *,
    strategies: Sequence[StrategyProposal],
    condition_evidence_by_strategy: Mapping[str, Any],
    selectivity_assessment_by_strategy: Mapping[
        str, ConditionSelectivityRepairAssessment | None
    ] | None = None,
) -> tuple[SingleStepRepairProposal, ...]:
    """Enumerate only existing realizations, strategies, and recipes."""

    assessments = selectivity_assessment_by_strategy or {}
    source_condition = condition_evidence_by_strategy.get(source_strategy.strategy_id)
    source_assessment = assessments.get(source_strategy.strategy_id)
    source_issues = collect_single_step_refinement_issues(
        source_strategy,
        condition_evidence=source_condition,
        condition_selectivity_assessment=source_assessment,
    )
    objective = _OBJECTIVE_BY_ISSUE.get(issue.kind, "find_alternative_strategy")
    source_realization_id = (
        source_strategy.representative.realization_id or source_strategy.strategy_id
    )

    def build(
        kind: SingleStepRepairKind,
        status: SingleStepRepairStatus,
        reason: str,
        *,
        target_strategy_id: str | None = None,
        target_realization_id: str | None = None,
        assessment_id: str | None = None,
        details: tuple[tuple[str, Any], ...] = (),
    ) -> SingleStepRepairProposal:
        identity = {
            "issue_id": issue.issue_id,
            "source_strategy_id": source_strategy.strategy_id,
            "source_realization_id": source_realization_id,
            "kind": kind,
            "status": status,
            "target_strategy_id": target_strategy_id,
            "target_realization_id": target_realization_id,
            "assessment_id": assessment_id,
            "schema_version": SINGLE_STEP_REFINEMENT_SCHEMA_VERSION,
        }
        return SingleStepRepairProposal(
            proposal_id=_stable_id("SSPROP1", identity),
            issue_id=issue.issue_id,
            source_strategy_id=source_strategy.strategy_id,
            source_realization_id=source_realization_id,
            repair_kind=kind,
            status=status,
            objective=objective,
            reason=reason,
            target_strategy_id=target_strategy_id,
            target_realization_id=target_realization_id,
            assessment_id=assessment_id,
            details=details,
        )

    proposals: list[SingleStepRepairProposal] = []
    if issue.kind == "selectivity":
        subject = _synthetic_route(
            source_strategy,
            condition_evidence=source_condition,
            condition_selectivity_assessment=source_assessment,
        ).steps[0]
        for assessment in assess_condition_selectivity_repairs(subject):
            proposals.append(
                build(
                    "condition_selectivity",
                    "actionable" if assessment.actionable else "unavailable",
                    assessment.reason,
                    assessment_id=(assessment.assessment_id if assessment.actionable else None),
                    details=(
                        ("assessment_id", assessment.assessment_id),
                        ("recipe_id", assessment.recipe_id),
                        ("recipe_core_id", assessment.recipe_core_id),
                        ("assessment_status", assessment.status),
                        (
                            "exact_condition_reference_ids",
                            assessment.exact_condition_reference_ids,
                        ),
                    ),
                )
            )

    for candidate in source_strategy.alternate_realizations:
        condition = (
            source_condition
            if _reaction_identity(candidate)
            == _reaction_identity(source_strategy.representative)
            else None
        )
        target_issues = collect_single_step_refinement_issues(
            source_strategy,
            candidate=candidate,
            condition_evidence=condition,
        )
        improved = _is_improvement(source_issues, target_issues, issue.kind)
        proposals.append(
            build(
                "alternate_realization",
                "actionable" if improved else "unavailable",
                (
                    "An already validated realization reduces the selected issue "
                    "without increasing strong or total deterministic issues."
                    if improved
                    else "The retained alternate realization does not verify as an issue-count improvement."
                ),
                target_strategy_id=(source_strategy.strategy_id if improved else None),
                target_realization_id=(candidate.realization_id if improved else None),
                details=(("candidate_realization_id", candidate.realization_id),),
            )
        )

    for target in strategies:
        if target.strategy_id == source_strategy.strategy_id:
            continue
        target_issues = collect_single_step_refinement_issues(
            target,
            condition_evidence=condition_evidence_by_strategy.get(target.strategy_id),
            condition_selectivity_assessment=assessments.get(target.strategy_id),
        )
        improved = _is_improvement(source_issues, target_issues, issue.kind)
        proposals.append(
            build(
                "alternate_strategy",
                "actionable" if improved else "unavailable",
                (
                    "An already validated strategy reduces the selected issue without "
                    "increasing strong or total deterministic issues."
                    if improved
                    else "The retained alternative strategy does not verify as an issue-count improvement."
                ),
                target_strategy_id=(target.strategy_id if improved else None),
                target_realization_id=(
                    target.representative.realization_id if improved else None
                ),
                details=(("candidate_strategy_id", target.strategy_id),),
            )
        )
    return tuple(proposals)


@dataclass(frozen=True)
class SingleStepVerificationGate:
    gate: str
    status: Literal["pass", "caution", "fail"]
    message: str


@dataclass(frozen=True)
class SingleStepVerificationReport:
    strategy_id: str
    realization_id: str
    status: Literal["verified", "verified_with_cautions", "failed"]
    gates: tuple[SingleStepVerificationGate, ...]
    issue_count: int
    strong_issue_count: int
    schema_version: str = SINGLE_STEP_REFINEMENT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def verify_single_step_strategy(
    strategy: StrategyProposal,
    *,
    candidate: GenericDisconnectionCandidate | None = None,
    condition_evidence: Any = None,
    condition_selectivity_assessment: ConditionSelectivityRepairAssessment
    | None = None,
) -> SingleStepVerificationReport:
    """Independently recheck identities and deterministic issue gates."""

    selected = candidate or strategy.representative
    gates: list[SingleStepVerificationGate] = []
    gates.append(
        SingleStepVerificationGate(
            "forward_validation",
            "pass" if selected.forward_validation_status == "verified_signature" else "fail",
            f"forward validation status is {selected.forward_validation_status}",
        )
    )
    split = split_reaction_smiles(selected.proposed_reaction_smiles)
    if split is None:
        gates.append(
            SingleStepVerificationGate(
                "reaction_parse", "fail", "proposed reaction SMILES could not be parsed"
            )
        )
    else:
        reactants, product = split
        canonical_product = canonical_smiles(product)
        canonical_target = canonical_smiles(selected.target_smiles)
        canonical_reactants = canonical_smiles(reactants)
        canonical_precursors = canonical_smiles(selected.precursor_smiles)
        gates.extend(
            (
                SingleStepVerificationGate(
                    "reaction_parse", "pass", "proposed reaction SMILES parsed"
                ),
                SingleStepVerificationGate(
                    "target_identity",
                    "pass" if canonical_product == canonical_target and canonical_target else "fail",
                    "reaction product matches the canonical target"
                    if canonical_product == canonical_target and canonical_target
                    else "reaction product does not match the canonical target",
                ),
                SingleStepVerificationGate(
                    "precursor_identity",
                    "pass"
                    if canonical_reactants == canonical_precursors and canonical_precursors
                    else "fail",
                    "reaction reactants match the canonical precursor set"
                    if canonical_reactants == canonical_precursors and canonical_precursors
                    else "reaction reactants do not match the canonical precursor set",
                ),
            )
        )
    issues = collect_single_step_refinement_issues(
        strategy,
        candidate=selected,
        condition_evidence=condition_evidence,
        condition_selectivity_assessment=condition_selectivity_assessment,
    )
    strong = sum(item.severity == "strong" for item in issues)
    gates.append(
        SingleStepVerificationGate(
            "deterministic_issue_profile",
            "caution" if issues else "pass",
            f"{len(issues)} unresolved deterministic issue(s); {strong} strong",
        )
    )
    if any(gate.status == "fail" for gate in gates):
        status = "failed"
    elif any(gate.status == "caution" for gate in gates):
        status = "verified_with_cautions"
    else:
        status = "verified"
    return SingleStepVerificationReport(
        strategy_id=strategy.strategy_id,
        realization_id=selected.realization_id or strategy.strategy_id,
        status=status,
        gates=tuple(gates),
        issue_count=len(issues),
        strong_issue_count=strong,
    )


__all__ = [
    "SINGLE_STEP_REFINEMENT_SCHEMA_VERSION",
    "SingleStepRepairProposal",
    "SingleStepVerificationGate",
    "SingleStepVerificationReport",
    "collect_single_step_refinement_issues",
    "enumerate_single_step_repair_proposals",
    "verify_single_step_strategy",
]
