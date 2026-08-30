"""Evidence-gated workbench for provider-backed real-world route examples."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Optional

from .generic_models import GenericTemplateLibrary
from .multistep import (
    ConditionEvidenceEvaluator,
    LiteratureLookup,
    MultistepRetrosynthesisResult,
    MultistepRetrosynthesisRoute,
    plan_multistep_routes,
)
from .provider_multistep import (
    ProviderBackedOneStepExpander,
    expansion_state_from_route,
)
from .route_refinement import (
    RouteRefinementIssue,
    RouteRepairProposal,
    collect_route_refinement_issues,
    enumerate_route_repair_proposals,
)
from .route_verification import RouteVerificationReport, verify_planned_route
from .step_precedents import StepPrecedentLookupResult, lookup_step_precedents
from .transition_orchestration import (
    ProviderExpansionOutcome,
    TransitionProviderOrchestrator,
)


ROUTE_WORKBENCH_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class RouteWorkbenchSettings:
    """Bounded deterministic search settings for one workbench run."""

    max_depth: int = 3
    molecular_weight_threshold: float = 150.0
    top_k_routes: int = 3
    per_step_top_k: int = 5
    beam_width: int = 12
    max_expansions: int = 20

    def __post_init__(self) -> None:
        if self.max_depth not in {2, 3}:
            raise ValueError("workbench route depth must be two or three")
        if self.molecular_weight_threshold <= 0:
            raise ValueError("workbench molecular-weight threshold must be positive")
        for value, name in (
            (self.top_k_routes, "top-k routes"),
            (self.per_step_top_k, "per-step top-k"),
            (self.beam_width, "beam width"),
            (self.max_expansions, "maximum expansions"),
        ):
            if value < 1:
                raise ValueError(f"{name} must be positive")


@dataclass(frozen=True)
class RouteWorkbenchStepEvidence:
    """Provider, precedent, and issue evidence for one retained route step."""

    step_id: str
    depth: int
    product_smiles: str
    precursor_smiles: tuple[str, ...]
    provider_id: str
    provider_rank: Optional[int]
    transition_id: Optional[str]
    operator_id: str
    strategy_id: str
    strategic_class: str
    strategic_complexity_score: float
    precursor_compatibility_disposition: str
    reaction_compatibility_disposition: str
    selectivity_warning_count: int
    condition_status: Optional[str]
    precedent_evidence: Optional[StepPrecedentLookupResult]
    issue_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible step evidence."""

        value = asdict(self)
        value["precedent_evidence"] = (
            self.precedent_evidence.to_dict()
            if self.precedent_evidence is not None
            else None
        )
        return value


@dataclass(frozen=True)
class RouteWorkbenchRouteReport:
    """One route plus deterministic verification and weakest-link evidence."""

    route: MultistepRetrosynthesisRoute
    verification: RouteVerificationReport
    step_evidence: tuple[RouteWorkbenchStepEvidence, ...]
    issues: tuple[RouteRefinementIssue, ...]
    weakest_issue_id: Optional[str]
    repair_proposals: tuple[RouteRepairProposal, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible route workbench evidence."""

        expansion_state = expansion_state_from_route(self.route)
        return {
            "route": self.route.to_dict(),
            "verification": self.verification.to_dict(),
            "expansion_state": expansion_state.to_dict(),
            "step_evidence": [item.to_dict() for item in self.step_evidence],
            "issues": [item.to_dict() for item in self.issues],
            "weakest_issue_id": self.weakest_issue_id,
            "repair_proposals": [item.to_dict() for item in self.repair_proposals],
        }


@dataclass(frozen=True)
class ProviderRouteWorkbenchResult:
    """Search result and evidence packets for one provider-backed target."""

    target_smiles: str
    provider_id: str
    settings: RouteWorkbenchSettings
    search_result: MultistepRetrosynthesisResult
    route_kind: str
    routes: tuple[RouteWorkbenchRouteReport, ...]
    provider_outcomes: tuple[ProviderExpansionOutcome, ...]
    schema_version: str = ROUTE_WORKBENCH_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a complete JSON-compatible workbench report."""

        return {
            "target_smiles": self.target_smiles,
            "provider_id": self.provider_id,
            "settings": asdict(self.settings),
            "search_result": self.search_result.to_dict(),
            "route_kind": self.route_kind,
            "routes": [item.to_dict() for item in self.routes],
            "provider_outcomes": [
                item.to_dict() for item in self.provider_outcomes
            ],
            "schema_version": self.schema_version,
        }


_ISSUE_KIND_PRIORITY = {
    "reaction_compatibility": 0,
    "precursor_compatibility": 1,
    "unresolved_leaf": 2,
    "selectivity": 3,
    "condition_gap": 4,
    "tactical_step": 5,
}


def _weakest_issue(
    issues: tuple[RouteRefinementIssue, ...],
) -> Optional[RouteRefinementIssue]:
    if not issues:
        return None
    return min(
        issues,
        key=lambda item: (
            0 if item.severity == "strong" else 1,
            _ISSUE_KIND_PRIORITY[item.kind],
            item.step_index if item.step_index is not None else 10_000,
            item.issue_id,
        ),
    )


def _step_precedents(
    route: MultistepRetrosynthesisRoute,
    library: GenericTemplateLibrary,
    step_index: int,
) -> Optional[StepPrecedentLookupResult]:
    try:
        return lookup_step_precedents(route.steps[step_index], library)
    except ValueError:
        return None


def _route_report(
    route: MultistepRetrosynthesisRoute,
    library: GenericTemplateLibrary,
    expander: ProviderBackedOneStepExpander,
) -> RouteWorkbenchRouteReport:
    issues = collect_route_refinement_issues(route)
    issue_ids_by_step: dict[str, list[str]] = {}
    for issue in issues:
        if issue.subject_type == "step":
            issue_ids_by_step.setdefault(issue.subject_id, []).append(issue.issue_id)
    step_evidence = []
    for index, step in enumerate(route.steps):
        attribution = expander.attribution_for_candidate(step.candidate)
        step_evidence.append(
            RouteWorkbenchStepEvidence(
                step_id=step.step_id,
                depth=step.depth,
                product_smiles=step.product_smiles,
                precursor_smiles=step.precursor_smiles,
                provider_id=expander.provider_id,
                provider_rank=(
                    attribution.provider_rank if attribution is not None else None
                ),
                transition_id=(
                    attribution.transition_id if attribution is not None else None
                ),
                operator_id=step.candidate.operator_id,
                strategy_id=step.candidate.strategy_id,
                strategic_class=step.candidate.strategic_class,
                strategic_complexity_score=(
                    step.candidate.strategic_complexity_score
                ),
                precursor_compatibility_disposition=(
                    step.candidate.precursor_compatibility_disposition
                ),
                reaction_compatibility_disposition=(
                    step.candidate.reaction_compatibility_disposition
                ),
                selectivity_warning_count=len(step.candidate.selectivity_warnings),
                condition_status=(
                    step.condition_evidence.status
                    if step.condition_evidence is not None
                    else None
                ),
                precedent_evidence=_step_precedents(route, library, index),
                issue_ids=tuple(sorted(issue_ids_by_step.get(step.step_id, ()))),
            )
        )
    weakest = _weakest_issue(issues)
    repairs = (
        enumerate_route_repair_proposals(route, weakest)
        if weakest is not None
        else ()
    )
    return RouteWorkbenchRouteReport(
        route=route,
        verification=verify_planned_route(route),
        step_evidence=tuple(step_evidence),
        issues=issues,
        weakest_issue_id=weakest.issue_id if weakest is not None else None,
        repair_proposals=repairs,
    )


def run_provider_route_workbench(
    target_smiles: str,
    library: GenericTemplateLibrary,
    literature_index: LiteratureLookup,
    orchestrator: TransitionProviderOrchestrator,
    provider_id: str,
    *,
    settings: RouteWorkbenchSettings | None = None,
    condition_evidence_evaluator: ConditionEvidenceEvaluator | None = None,
) -> ProviderRouteWorkbenchResult:
    """Run bounded multistep search and compose canonical evidence tools."""

    active = settings or RouteWorkbenchSettings()
    expander = ProviderBackedOneStepExpander(orchestrator, provider_id)
    search = plan_multistep_routes(
        target_smiles,
        library,
        literature_index,
        max_depth=active.max_depth,
        molecular_weight_threshold=active.molecular_weight_threshold,
        top_k_routes=active.top_k_routes,
        per_step_top_k=active.per_step_top_k,
        beam_width=active.beam_width,
        max_expansions=active.max_expansions,
        condition_evidence_evaluator=condition_evidence_evaluator,
        expander=expander,
    )
    retained = search.routes if search.routes else search.partial_routes
    route_kind = "solved" if search.routes else "partial"
    reports = tuple(
        _route_report(route, library, expander) for route in retained
    )
    return ProviderRouteWorkbenchResult(
        target_smiles=search.target_smiles,
        provider_id=provider_id,
        settings=active,
        search_result=search,
        route_kind=route_kind,
        routes=reports,
        provider_outcomes=expander.outcomes,
    )


__all__ = [
    "ROUTE_WORKBENCH_SCHEMA_VERSION",
    "ProviderRouteWorkbenchResult",
    "RouteWorkbenchRouteReport",
    "RouteWorkbenchSettings",
    "RouteWorkbenchStepEvidence",
    "run_provider_route_workbench",
]
