"""Closed application capabilities over authoritative domain services."""

from __future__ import annotations

from dataclasses import dataclass, replace
from types import SimpleNamespace
from typing import Any, Dict, Mapping, Tuple

from condition_recommender.models import GenericRecommendationResult
from core_retrosynthesis import (
    MultistepRetrosynthesisRoute,
    RouteCandidateExclusion,
    RouteRepairProposal,
    RouteRefinementIntent,
    RouteSearchPolicyDelta,
    apply_route_search_delta,
    assess_condition_selectivity_repairs,
    build_route_refinement_plan,
    collect_route_refinement_issues,
    collect_single_step_refinement_issues,
    enumerate_single_step_repair_proposals,
    enumerate_route_repair_proposals,
    lookup_step_precedents,
    summarize_route_refinement,
    verify_planned_route,
    verify_single_step_strategy,
)
from reactive_taxonomy import audit_target

from chem_coworker.contracts import (
    ConditionRequest,
    ConditionReviewSettings,
    MultistepRetrosynthesisRequest,
    MultistepRetrosynthesisResponse,
    RetrosynthesisRequest,
    RetrosynthesisResponse,
)
from chem_coworker.multistep import MultistepRetrosynthesisCoworker
from chem_coworker.retrosynthesis import RetrosynthesisCoworker
from chem_coworker.retrosynthesis_rendering import render_retrosynthesis
from chem_coworker.service import ConditionCoworker

from .contracts import (
    AssistanceRequest,
    EvidenceItem,
    canonical_json,
    stable_assistance_id,
)
from .evidence import (
    ConditionEvidenceProjection,
    MultistepEvidenceProjection,
    RetrosynthesisEvidenceProjection,
    project_condition_result,
    project_multistep_response,
    project_retrosynthesis_response,
)


@dataclass(frozen=True)
class CapabilityResult:
    """Validated evidence added by one capability invocation."""

    result_ref: str
    evidence: Tuple[EvidenceItem, ...]
    packet: Mapping[str, Any]
    authoritative_result: object | None = None
    register_result_ref: bool = True


class ChemistryCapabilities:
    """Read-only molecular tools shared by retrosynthesis assistance modes."""

    def audit_target(self, request: AssistanceRequest) -> CapabilityResult:
        """Return compact RDKit-backed target evidence without changing it."""

        if request.mode not in {"retro", "multistep"}:
            raise ValueError("target audit requires a retrosynthesis mode")
        audit = audit_target(request.structure_input)
        result_ref = stable_assistance_id("TARGETAUDIT", audit.to_dict())
        evidence = EvidenceItem(
            evidence_id="target.audit",
            layer="observation",
            source_id=result_ref,
            payload_type="target_audit",
            payload=audit.to_dict(),
            provenance="deterministic_inference",
            schema_versions={"target_audit": audit.schema_version},
            uncertainty=(
                audit.error
                or "; ".join(audit.warnings)
                or "Reactive sites are graph-derived hypotheses, not reaction claims."
            ),
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=(evidence,),
            packet={
                "result_ref": result_ref,
                "evidence": [_evidence_packet(evidence)],
            },
            authoritative_result=audit,
            register_result_ref=False,
        )


class ConditionCapabilities:
    """Condition assistance adapters with request-local result storage."""

    def __init__(self, coworker: ConditionCoworker) -> None:
        self._coworker = coworker
        self._results: Dict[str, GenericRecommendationResult] = {}
        self._projections: Dict[str, ConditionEvidenceProjection] = {}

    def recommend_conditions(self, request: AssistanceRequest) -> CapabilityResult:
        """Run the canonical recommender using only confirmed constraints."""

        if request.mode != "conditions":
            raise ValueError("condition recommendation requires conditions mode")
        kwargs: Dict[str, Any] = {
            "reaction_smiles": request.structure_input,
            "condition_constraints": request.condition_constraints,
            "review": ConditionReviewSettings(mode="off"),
        }
        for constraint in request.constraints:
            if constraint.owner != "condition_recommender":
                continue
            if constraint.provenance == "model_proposed":
                continue
            if constraint.kind in {
                "top_k",
                "minimum_pool_size",
                "unrestricted_fallback",
                "ranking_profile",
            }:
                kwargs[constraint.kind] = constraint.value
            elif constraint.kind == "ranking_weights":
                kwargs[constraint.kind] = dict(constraint.value)
        response = self._coworker.recommend(ConditionRequest(**kwargs))
        projection = project_condition_result(response.result)
        self._results[projection.result_ref] = response.result
        self._projections[projection.result_ref] = projection
        initial_ids = {
            item["evidence_id"] for item in projection.initial_packet()["evidence"]
        }
        evidence = tuple(
            item for item in projection.evidence if item.evidence_id in initial_ids
        )
        return CapabilityResult(
            result_ref=projection.result_ref,
            evidence=evidence,
            packet=projection.initial_packet(),
            authoritative_result=response.result,
        )

    def inspect_condition_candidate(
        self,
        *,
        result_ref: str,
        candidate_alias: str,
    ) -> CapabilityResult:
        """Expose full existing evidence for one known recommendation."""

        projection = self._projection(result_ref)
        evidence = tuple(
            item
            for item in projection.candidate_evidence(candidate_alias)
            if not item.evidence_id.endswith(".summary")
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=evidence,
            packet=projection.inspection_packet(candidate_alias),
            authoritative_result=self._results[result_ref],
        )

    def compare_condition_candidates(
        self,
        *,
        result_ref: str,
        candidate_aliases: Tuple[str, ...],
    ) -> CapabilityResult:
        """Expose existing candidate differences without reranking chemistry."""

        projection = self._projection(result_ref)
        packet = projection.comparison_packet(candidate_aliases)
        evidence_ids = {
            item["evidence_id"] for item in packet["evidence"]
        }
        return CapabilityResult(
            result_ref=result_ref,
            evidence=tuple(
                item for item in projection.evidence if item.evidence_id in evidence_ids
            ),
            packet=packet,
            authoritative_result=self._results[result_ref],
        )

    def result(self, result_ref: str) -> GenericRecommendationResult:
        """Return an authoritative result by its application-managed reference."""

        try:
            return self._results[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale condition result: {result_ref}") from exc

    def request_fingerprint(self, request: AssistanceRequest) -> str:
        """Return the normalized domain input used for repeat detection."""

        return canonical_json(
            {
                "structure_input": request.structure_input,
                "constraints": request.constraints,
                "condition_constraints": request.condition_constraints,
            }
        )

    def _projection(self, result_ref: str) -> ConditionEvidenceProjection:
        try:
            return self._projections[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale condition result: {result_ref}") from exc


class RetrosynthesisCapabilities:
    """One-step adapters that expose only already validated strategies."""

    def __init__(self, coworker: RetrosynthesisCoworker) -> None:
        self._coworker = coworker
        self._responses: Dict[str, RetrosynthesisResponse] = {}
        self._projections: Dict[str, RetrosynthesisEvidenceProjection] = {}
        self._accepted_targets: Dict[str, frozenset[str]] = {}

    def disconnect_target(self, request: AssistanceRequest) -> CapabilityResult:
        if request.mode != "retro":
            raise ValueError("target disconnection requires retro mode")
        response = self._coworker.disconnect(
            RetrosynthesisRequest(
                target_smiles=request.structure_input,
                review=ConditionReviewSettings(mode="off"),
            )
        )
        projection = project_retrosynthesis_response(response)
        self._responses[projection.result_ref] = response
        self._projections[projection.result_ref] = projection
        self._accepted_targets[projection.result_ref] = frozenset()
        initial_ids = {
            item["evidence_id"] for item in projection.initial_packet()["evidence"]
        }
        return CapabilityResult(
            result_ref=projection.result_ref,
            evidence=tuple(
                item for item in projection.evidence if item.evidence_id in initial_ids
            ),
            packet=projection.initial_packet(),
            authoritative_result=response,
        )

    def inspect_strategy(self, result_ref: str, alias: str) -> CapabilityResult:
        projection = self._projection(result_ref)
        evidence = tuple(
            item
            for item in projection.strategy_evidence(alias)
            if not item.evidence_id.endswith(".summary")
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=evidence,
            packet=projection.inspection_packet(alias),
            authoritative_result=self._responses[result_ref],
        )

    def inspect_strategy_conditions(
        self,
        result_ref: str,
        alias: str,
    ) -> CapabilityResult:
        projection = self._projection(result_ref)
        evidence = tuple(
            item
            for item in projection.strategy_evidence(alias)
            if item.evidence_id.endswith(".conditions")
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=evidence,
            packet=_packet(result_ref, evidence, strategy_alias=alias),
            authoritative_result=self._responses[result_ref],
        )

    def compare_strategies(
        self,
        result_ref: str,
        aliases: Tuple[str, ...],
    ) -> CapabilityResult:
        projection = self._projection(result_ref)
        packet = projection.comparison_packet(aliases)
        evidence_ids = {item["evidence_id"] for item in packet["evidence"]}
        return CapabilityResult(
            result_ref=result_ref,
            evidence=tuple(
                item for item in projection.evidence if item.evidence_id in evidence_ids
            ),
            packet=packet,
            authoritative_result=self._responses[result_ref],
        )

    def apply_repair(
        self,
        request: AssistanceRequest,
        result_ref: str,
        *,
        proposal_id: str,
    ) -> CapabilityResult:
        """Apply one ID-selected recipe, realization, or retained strategy."""

        response = self.result(result_ref)
        projection = self._projection(result_ref)
        proposal_evidence = next(
            (
                item
                for item in projection.evidence
                if item.payload_type == "single_step_repair_proposal"
                and item.payload.get("proposal_id") == proposal_id
            ),
            None,
        )
        if proposal_evidence is None:
            raise ValueError("unknown deterministic single-step repair proposal")
        if proposal_evidence.payload.get("status") != "actionable":
            raise ValueError("single-step repair proposal is not actionable")
        source_id = str(proposal_evidence.payload.get("source_strategy_id") or "")
        source = next(
            (item for item in response.strategies if item.strategy_id == source_id),
            None,
        )
        if source is None:
            raise ValueError("single-step repair source strategy is unavailable")
        conditions = {item.strategy_id: item for item in response.condition_evidence}
        source_condition = conditions.get(source.strategy_id)
        source_issues = collect_single_step_refinement_issues(
            source,
            condition_evidence=(source_condition.evidence if source_condition else None),
            condition_selectivity_assessment=(
                source_condition.condition_selectivity_assessment
                if source_condition
                else None
            ),
        )
        issue_id = str(proposal_evidence.payload.get("issue_id") or "")
        issue = next((item for item in source_issues if item.issue_id == issue_id), None)
        if issue is None:
            raise ValueError("single-step repair issue is no longer present")
        proposals = enumerate_single_step_repair_proposals(
            source,
            issue,
            strategies=response.strategies,
            condition_evidence_by_strategy={
                key: value.evidence for key, value in conditions.items()
            },
            selectivity_assessment_by_strategy={
                key: value.condition_selectivity_assessment
                for key, value in conditions.items()
            },
        )
        proposal = next(
            (
                item
                for item in proposals
                if item.proposal_id == proposal_id and item.status == "actionable"
            ),
            None,
        )
        if proposal is None:
            raise ValueError("single-step repair proposal is no longer actionable")
        target_key = (
            proposal.assessment_id
            or proposal.target_realization_id
            or proposal.target_strategy_id
            or proposal.proposal_id
        )
        if target_key in self._accepted_targets.get(result_ref, frozenset()):
            raise ValueError("single-step repair would revisit an accepted target")

        updated_response = response
        target_strategy = source
        target_condition = source_condition
        if proposal.repair_kind == "condition_selectivity":
            subject = SimpleNamespace(
                step_id=(
                    source.representative.realization_id or source.strategy_id
                ),
                candidate=source.representative,
                condition_evidence=(
                    source_condition.evidence if source_condition else None
                ),
                condition_selectivity_assessment=(
                    source_condition.condition_selectivity_assessment
                    if source_condition
                    else None
                ),
            )
            assessment = next(
                (
                    item
                    for item in assess_condition_selectivity_repairs(subject)
                    if item.assessment_id == proposal.assessment_id and item.actionable
                ),
                None,
            )
            if assessment is None or source_condition is None:
                raise ValueError("condition-selectivity evidence is no longer actionable")
            target_condition = replace(
                source_condition,
                condition_selectivity_assessment=assessment,
            )
            conditions[source.strategy_id] = target_condition
            updated_response = replace(
                response,
                condition_evidence=tuple(
                    conditions.get(item.strategy_id, item)
                    for item in response.condition_evidence
                ),
                selected_strategy_id=source.strategy_id,
                selected_realization_id=(
                    source.representative.realization_id or source.strategy_id
                ),
            )
        elif proposal.repair_kind == "alternate_strategy":
            target_strategy = next(
                (
                    item
                    for item in response.strategies
                    if item.strategy_id == proposal.target_strategy_id
                ),
                None,
            )
            if target_strategy is None:
                raise ValueError("target single-step strategy is unavailable")
            target_condition = conditions.get(target_strategy.strategy_id)
            updated_response = replace(
                response,
                selected_strategy_id=target_strategy.strategy_id,
                selected_realization_id=(
                    target_strategy.representative.realization_id
                    or target_strategy.strategy_id
                ),
            )
        else:
            target_candidate = next(
                (
                    item
                    for item in source.alternate_realizations
                    if item.realization_id == proposal.target_realization_id
                ),
                None,
            )
            if target_candidate is None:
                raise ValueError("target single-step realization is unavailable")
            retained_alternates = (
                source.representative,
                *(
                    item
                    for item in source.alternate_realizations
                    if item.realization_id != target_candidate.realization_id
                ),
            )
            target_strategy = replace(
                source,
                representative=target_candidate,
                alternate_realizations=retained_alternates,
            )
            replacement_conditions = self._coworker.recommend_strategy_conditions(
                target_strategy,
                response.request,
            )
            target_condition = replacement_conditions
            strategies = tuple(
                target_strategy if item.strategy_id == source.strategy_id else item
                for item in response.strategies
            )
            updated_conditions = tuple(
                item
                for item in response.condition_evidence
                if item.strategy_id != source.strategy_id
            )
            if replacement_conditions is not None:
                updated_conditions += (replacement_conditions,)
            updated_response = replace(
                response,
                strategies=strategies,
                condition_evidence=updated_conditions,
                selected_strategy_id=target_strategy.strategy_id,
                selected_realization_id=(
                    target_candidate.realization_id or target_strategy.strategy_id
                ),
            )

        target_issues = collect_single_step_refinement_issues(
            target_strategy,
            condition_evidence=(target_condition.evidence if target_condition else None),
            condition_selectivity_assessment=(
                target_condition.condition_selectivity_assessment
                if target_condition
                else None
            ),
        )
        source_relevant = sum(item.kind == issue.kind for item in source_issues)
        target_relevant = sum(item.kind == issue.kind for item in target_issues)
        source_strong = sum(item.severity == "strong" for item in source_issues)
        target_strong = sum(item.severity == "strong" for item in target_issues)
        verification = verify_single_step_strategy(
            target_strategy,
            condition_evidence=(target_condition.evidence if target_condition else None),
            condition_selectivity_assessment=(
                target_condition.condition_selectivity_assessment
                if target_condition
                else None
            ),
        )
        accepted = (
            verification.status != "failed"
            and target_relevant < source_relevant
            and target_strong <= source_strong
            and len(target_issues) <= len(source_issues)
        )
        retained_response = updated_response if accepted else response
        if accepted:
            retained_response = replace(
                retained_response,
                answer=render_retrosynthesis(retained_response),
            )
        retained_projection = project_retrosynthesis_response(retained_response)
        outcome_id = stable_assistance_id(
            "SSREFINE",
            {
                "source_result_ref": result_ref,
                "proposal_id": proposal.proposal_id,
                "accepted": accepted,
                "retained_result_ref": retained_projection.result_ref,
            },
        )
        outcome = EvidenceItem(
            evidence_id=f"refinement.{outcome_id}",
            layer="route",
            source_id=result_ref,
            payload_type="single_step_refinement_outcome",
            payload={
                "proposal_id": proposal.proposal_id,
                "repair_kind": proposal.repair_kind,
                "status": (
                    "improved_alternative_found"
                    if accepted
                    else "alternatives_found_no_verified_improvement"
                ),
                "source_strategy_id": source.strategy_id,
                "retained_strategy_id": (
                    target_strategy.strategy_id if accepted else source.strategy_id
                ),
                "source_result_ref": result_ref,
                "retained_result_ref": retained_projection.result_ref,
                "source_relevant_issue_count": source_relevant,
                "retained_relevant_issue_count": (
                    target_relevant if accepted else source_relevant
                ),
                "source_issue_count": len(source_issues),
                "retained_issue_count": (
                    len(target_issues) if accepted else len(source_issues)
                ),
                "source_preserved": True,
                "schema_version": proposal.schema_version,
            },
            provenance="deterministic_inference",
            schema_versions={"single_step_refinement": proposal.schema_version},
            uncertainty=(
                "Deterministic issue reduction is not experimental proof of yield."
            ),
        )
        verification_evidence = EvidenceItem(
            evidence_id=f"refinement.{outcome_id}.verification",
            layer="route",
            source_id=target_strategy.strategy_id,
            payload_type="single_step_refinement_verification",
            payload={"accepted": accepted, **verification.to_dict()},
            provenance="deterministic_inference",
            schema_versions={
                "single_step_refinement": verification.schema_version
            },
            uncertainty=(
                "Deterministic verification is not experimental proof of success."
            ),
        )
        if accepted:
            history = frozenset(
                {*self._accepted_targets.get(result_ref, frozenset()), target_key}
            )
            self._responses[retained_projection.result_ref] = retained_response
            self._projections[retained_projection.result_ref] = retained_projection
            self._accepted_targets[retained_projection.result_ref] = history
            initial_ids = {
                item["evidence_id"]
                for item in retained_projection.initial_packet()["evidence"]
            }
            exposed = tuple(
                item
                for item in retained_projection.evidence
                if item.evidence_id in initial_ids
            )
        else:
            exposed = ()
        return CapabilityResult(
            result_ref=retained_projection.result_ref,
            evidence=(outcome, verification_evidence, *exposed),
            packet={
                **retained_projection.initial_packet(),
                "refinement": outcome.payload,
                "verification": verification_evidence.payload,
            },
            authoritative_result=retained_response,
            register_result_ref=accepted,
        )

    def verify_strategy(self, result_ref: str, alias: str) -> CapabilityResult:
        """Independently verify one retained one-step strategy."""

        projection = self._projection(result_ref)
        strategy_id = dict(projection.strategy_aliases).get(alias)
        if strategy_id is None:
            raise ValueError(f"unknown retrosynthesis strategy alias: {alias}")
        response = self.result(result_ref)
        strategy = next(
            item for item in response.strategies if item.strategy_id == strategy_id
        )
        condition = next(
            (
                item
                for item in response.condition_evidence
                if item.strategy_id == strategy_id
            ),
            None,
        )
        report = verify_single_step_strategy(
            strategy,
            condition_evidence=(condition.evidence if condition else None),
            condition_selectivity_assessment=(
                condition.condition_selectivity_assessment if condition else None
            ),
        )
        evidence = EvidenceItem(
            evidence_id=f"{alias}.verification",
            layer="route",
            source_id=strategy_id,
            payload_type="single_step_verification",
            payload={"strategy_alias": alias, **report.to_dict()},
            provenance="deterministic_inference",
            schema_versions={"single_step_refinement": report.schema_version},
            uncertainty=(
                "Deterministic verification is not experimental proof of success."
            ),
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=(evidence,),
            packet={
                "result_ref": result_ref,
                "strategy_alias": alias,
                "evidence": [_evidence_packet(evidence)],
            },
            authoritative_result=response,
            register_result_ref=False,
        )

    def result(self, result_ref: str) -> RetrosynthesisResponse:
        try:
            return self._responses[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale retrosynthesis result: {result_ref}") from exc

    def _projection(self, result_ref: str) -> RetrosynthesisEvidenceProjection:
        try:
            return self._projections[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale retrosynthesis result: {result_ref}") from exc


class MultistepCapabilities:
    """Route adapters that preserve solved/partial status and search diagnostics."""

    def __init__(self, coworker: MultistepRetrosynthesisCoworker) -> None:
        self._coworker = coworker
        self._responses: Dict[str, MultistepRetrosynthesisResponse] = {}
        self._requests: Dict[str, MultistepRetrosynthesisRequest] = {}
        self._projections: Dict[str, MultistepEvidenceProjection] = {}
        self._candidate_exclusions: Dict[
            str, tuple[RouteCandidateExclusion, ...]
        ] = {}
        self._accepted_route_ids: Dict[str, frozenset[str]] = {}

    def plan_routes(self, request: AssistanceRequest) -> CapabilityResult:
        if request.mode != "multistep":
            raise ValueError("route planning requires multistep mode")
        policy = request.route_search_policy
        domain_request = MultistepRetrosynthesisRequest(
            target_smiles=request.structure_input,
            max_depth=policy.initial_max_depth,  # type: ignore[arg-type]
            beam_width=policy.initial_beam_width,
            max_expansions=policy.initial_max_expansions,
            review=ConditionReviewSettings(mode="off"),
        )
        return self._plan(domain_request)

    def retry_route_search(
        self,
        request: AssistanceRequest,
        result_ref: str,
        delta: RouteSearchPolicyDelta,
    ) -> CapabilityResult:
        prior = self._request(result_ref)
        depth, beam, expansions = apply_route_search_delta(
            request.route_search_policy,
            delta,
            current_max_depth=prior.max_depth,
            current_beam_width=prior.beam_width,
            current_max_expansions=prior.max_expansions,
        )
        expanded = MultistepRetrosynthesisRequest(
            **{
                **prior.__dict__,
                "max_depth": depth,
                "beam_width": beam,
                "max_expansions": expansions,
            }
        )
        return self._plan(
            expanded,
            candidate_exclusions=self._candidate_exclusions.get(result_ref, ()),
            accepted_route_ids=self._accepted_route_ids.get(result_ref, frozenset()),
        )

    def apply_repair(
        self,
        request: AssistanceRequest,
        result_ref: str,
        *,
        proposal_id: str,
    ) -> CapabilityResult:
        """Execute one actionable repair using only its deterministic identity."""

        projection = self._projection(result_ref)
        proposal_evidence = next(
            (
                item
                for item in projection.evidence
                if item.payload_type == "route_repair_proposal"
                and item.payload.get("proposal_id") == proposal_id
            ),
            None,
        )
        if proposal_evidence is None:
            raise ValueError("unknown deterministic repair proposal")
        if proposal_evidence.payload.get("status") != "actionable":
            raise ValueError("repair proposal is not actionable")
        route_alias = str(proposal_evidence.payload.get("route_alias") or "")
        step_index = int(proposal_evidence.payload.get("step_index") or 0)
        source_route = self._route(result_ref, route_alias)
        if step_index < 1 or step_index > len(source_route.steps):
            raise ValueError("repair proposal step index is out of range")
        issue_id = str(proposal_evidence.payload.get("issue_id") or "")
        source_issues = collect_route_refinement_issues(source_route)
        issue = next(
            (item for item in source_issues if item.issue_id == issue_id),
            None,
        )
        if issue is None:
            raise ValueError("repair proposal references an unavailable route issue")
        proposal = next(
            (
                item
                for item in enumerate_route_repair_proposals(source_route, issue)
                if item.proposal_id == proposal_id and item.status == "actionable"
            ),
            None,
        )
        if proposal is None or proposal.refinement_method is None:
            raise ValueError("repair proposal is no longer actionable")
        source_step = source_route.steps[step_index - 1]
        if proposal.refinement_method == "condition_selectivity":
            return self._apply_condition_selectivity_repair(
                result_ref,
                route_alias,
                source_route,
                source_step,
                proposal,
            )
        intent = RouteRefinementIntent(
            source_route_id=source_route.route_id,
            source_step_id=source_step.step_id,
            objective=proposal.objective,
            method=proposal.refinement_method,
            issue_ids=(issue.issue_id,),
            maximum_added_steps=proposal.maximum_added_steps,
        )
        plan = build_route_refinement_plan(source_route, intent, source_issues)
        prior = self._request(result_ref)
        policy = request.route_search_policy
        depth = min(
            policy.maximum_max_depth,
            prior.max_depth + proposal.maximum_added_steps,
        )
        refined_request = MultistepRetrosynthesisRequest(
            **{
                **prior.__dict__,
                "max_depth": depth,
                "beam_width": policy.maximum_beam_width,
                "max_expansions": policy.maximum_max_expansions,
            }
        )
        prior_exclusions = self._candidate_exclusions.get(result_ref, ())
        cumulative_exclusions = tuple(
            dict.fromkeys((*prior_exclusions, plan.exclusion))
        )
        response = self._coworker.plan(
            refined_request,
            candidate_exclusions=cumulative_exclusions,
        )
        refined_projection = project_multistep_response(response)
        if response.result is None:
            raise ValueError("deterministic refinement returned no route result")
        accepted_history = frozenset(
            {
                *self._accepted_route_ids.get(result_ref, frozenset()),
                source_route.route_id,
            }
        )
        outcome = summarize_route_refinement(
            source_route,
            intent,
            response.result,
            previously_accepted_route_ids=accepted_history,
        )
        accepted_route = next(
            (
                route
                for route in (*response.result.routes, *response.result.partial_routes)
                if route.route_id == outcome.accepted_route_id
            ),
            None,
        )
        accepted_alias = next(
            (
                alias
                for alias, route_id in refined_projection.route_aliases
                if route_id == outcome.accepted_route_id
            ),
            None,
        )
        retained_route = accepted_route or source_route
        retained_result_ref = (
            refined_projection.result_ref if accepted_route is not None else result_ref
        )
        retained_route_alias = accepted_alias or route_alias
        verification = verify_planned_route(retained_route)
        lineage = EvidenceItem(
            evidence_id=f"refinement.{intent.intent_id}",
            layer="route",
            source_id=refined_projection.result_ref,
            payload_type="route_refinement_outcome",
            payload={
                **outcome.to_dict(),
                "source_result_ref": result_ref,
                "refined_result_ref": refined_projection.result_ref,
                "retained_result_ref": retained_result_ref,
                "retained_route_alias": retained_route_alias,
                "excluded_candidate_scope": plan.excluded_candidate_scope,
                "cumulative_exclusion_count": len(cumulative_exclusions),
            },
            provenance="deterministic_inference",
            schema_versions={"route_refinement": outcome.schema_version},
            uncertainty=(
                "A lower deterministic issue count is not experimental proof of "
                "reaction success."
            ),
        )
        verification_evidence = EvidenceItem(
            evidence_id=f"refinement.{intent.intent_id}.verification",
            layer="route",
            source_id=retained_route.route_id,
            payload_type="route_refinement_verification",
            payload={
                "route_alias": retained_route_alias,
                "accepted": accepted_route is not None,
                **verification.to_dict(),
            },
            provenance="deterministic_inference",
            schema_versions={"route_verification": verification.schema_version},
            uncertainty=(
                "Deterministic route verification is not experimental proof of "
                "yield or operational feasibility."
            ),
        )
        next_history = frozenset(
            {
                *accepted_history,
                *(
                    (accepted_route.route_id,)
                    if accepted_route is not None
                    else ()
                ),
            }
        )
        self._store(
            refined_request,
            response,
            refined_projection,
            candidate_exclusions=cumulative_exclusions,
            accepted_route_ids=next_history,
        )
        # A failed attempt retains the authoritative source, but its exclusion is
        # still remembered so the next proposal cannot silently repeat the search.
        self._candidate_exclusions[result_ref] = cumulative_exclusions
        self._accepted_route_ids[result_ref] = accepted_history
        initial_ids = {
            item["evidence_id"]
            for item in refined_projection.initial_packet()["evidence"]
        }
        refined_summaries = tuple(
            item
            for item in refined_projection.evidence
            if item.evidence_id in initial_ids
        )
        exposed = refined_summaries if accepted_route is not None else ()
        retained_packet = (
            refined_projection.initial_packet()
            if accepted_route is not None
            else {
                "result_ref": result_ref,
                "route_aliases": [route_alias],
                "evidence": [],
            }
        )
        packet = {
            **retained_packet,
            "result_ref": retained_result_ref,
            "search_result_ref": refined_projection.result_ref,
            "refinement": lineage.payload,
            "verification": verification_evidence.payload,
        }
        return CapabilityResult(
            result_ref=retained_result_ref,
            evidence=(lineage, verification_evidence, *exposed),
            packet=packet,
            authoritative_result=(
                response if accepted_route is not None else self._responses[result_ref]
            ),
            register_result_ref=accepted_route is not None,
        )

    def _apply_condition_selectivity_repair(
        self,
        result_ref: str,
        route_alias: str,
        source_route: MultistepRetrosynthesisRoute,
        source_step: Any,
        proposal: RouteRepairProposal,
    ) -> CapabilityResult:
        """Attach one supported existing recipe without changing route chemistry."""

        details = dict(proposal.details)
        assessment_id = str(details.get("assessment_id") or "")
        assessment = next(
            (
                value
                for value in assess_condition_selectivity_repairs(source_step)
                if value.assessment_id == assessment_id and value.actionable
            ),
            None,
        )
        if assessment is None:
            raise ValueError(
                "condition-selectivity repair assessment is no longer actionable"
            )
        intent = RouteRefinementIntent(
            source_route_id=source_route.route_id,
            source_step_id=source_step.step_id,
            objective=proposal.objective,
            method="condition_selectivity",
            issue_ids=(proposal.issue_id,),
            maximum_added_steps=0,
        )
        updated_step = replace(
            source_step,
            condition_selectivity_assessment=assessment,
        )
        updated_steps = tuple(
            updated_step if step.step_id == source_step.step_id else step
            for step in source_route.steps
        )
        updated_route = replace(source_route, steps=updated_steps)
        source_response = self._responses[result_ref]
        if source_response.result is None:
            raise ValueError("multistep response has no route result")
        updated_result = replace(
            source_response.result,
            routes=tuple(
                updated_route if route.route_id == source_route.route_id else route
                for route in source_response.result.routes
            ),
            partial_routes=tuple(
                updated_route if route.route_id == source_route.route_id else route
                for route in source_response.result.partial_routes
            ),
        )
        updated_response = replace(source_response, result=updated_result)
        projection = project_multistep_response(updated_response)
        updated_alias = next(
            (
                alias
                for alias, route_id in projection.route_aliases
                if route_id == updated_route.route_id
            ),
            route_alias,
        )
        source_issue_count = len(collect_route_refinement_issues(source_route))
        updated_issue_count = len(collect_route_refinement_issues(updated_route))
        if updated_issue_count >= source_issue_count:
            raise ValueError(
                "condition-selectivity repair did not reduce deterministic issues"
            )
        verification = verify_planned_route(updated_route)
        hard_failures = tuple(
            gate.gate
            for gate in verification.gates
            if gate.gate
            in {
                "route_tree_integrity",
                "target_identity",
                "step_graph_consistency",
                "forward_validation",
            }
            and gate.status == "fail"
        )
        if hard_failures:
            raise ValueError(
                "condition-selectivity repair failed hard route verification"
            )
        lineage = EvidenceItem(
            evidence_id=f"refinement.{intent.intent_id}",
            layer="route",
            source_id=projection.result_ref,
            payload_type="route_refinement_outcome",
            payload={
                "schema_version": proposal.schema_version,
                "intent": intent.to_dict(),
                "status": "improved_alternative_found",
                "source_route_id": source_route.route_id,
                "accepted_route_id": updated_route.route_id,
                "retained_route_id": updated_route.route_id,
                "source_route_preserved": True,
                "source_result_ref": result_ref,
                "refined_result_ref": projection.result_ref,
                "retained_result_ref": projection.result_ref,
                "retained_route_alias": updated_alias,
                "repair_kind": proposal.repair_kind,
                "cumulative_exclusion_count": len(
                    self._candidate_exclusions.get(result_ref, ())
                ),
                "source_issue_count": source_issue_count,
                "retained_issue_count": updated_issue_count,
                "condition_selectivity_assessment": assessment.to_dict(),
            },
            provenance="deterministic_inference",
            schema_versions={"route_refinement": proposal.schema_version},
            uncertainty=(
                "Condition-supported endpoint preference is not experimental proof "
                "of selectivity or yield."
            ),
        )
        verification_evidence = EvidenceItem(
            evidence_id=f"refinement.{intent.intent_id}.verification",
            layer="route",
            source_id=updated_route.route_id,
            payload_type="route_refinement_verification",
            payload={
                "route_alias": updated_alias,
                "accepted": True,
                **verification.to_dict(),
            },
            provenance="deterministic_inference",
            schema_versions={"route_verification": verification.schema_version},
            uncertainty=(
                "Deterministic route verification is not experimental proof of "
                "yield or operational feasibility."
            ),
        )
        accepted_ids = frozenset(
            {
                *self._accepted_route_ids.get(result_ref, frozenset()),
                updated_route.route_id,
            }
        )
        self._store(
            self._request(result_ref),
            updated_response,
            projection,
            candidate_exclusions=self._candidate_exclusions.get(result_ref, ()),
            accepted_route_ids=accepted_ids,
        )
        initial_ids = {
            item["evidence_id"] for item in projection.initial_packet()["evidence"]
        }
        exposed = tuple(
            item for item in projection.evidence if item.evidence_id in initial_ids
        )
        return CapabilityResult(
            result_ref=projection.result_ref,
            evidence=(lineage, verification_evidence, *exposed),
            packet={
                **projection.initial_packet(),
                "refinement": lineage.payload,
                "verification": verification_evidence.payload,
            },
            authoritative_result=updated_response,
        )

    def search_step_precedents(
        self,
        result_ref: str,
        alias: str,
        *,
        step_index: int,
    ) -> CapabilityResult:
        """Expose admitted precedents supporting one existing route step."""

        route = self._route(result_ref, alias)
        if step_index < 1 or step_index > len(route.steps):
            raise ValueError("route precedent step index is out of range")
        result = lookup_step_precedents(
            route.steps[step_index - 1],
            self._coworker.library,
        )
        prefix = f"{alias}.step-{step_index}.precedents"
        evidence = [
            EvidenceItem(
                evidence_id=f"{prefix}.summary",
                layer="route",
                source_id=result_ref,
                payload_type="step_precedent_lookup",
                payload={
                    "route_alias": alias,
                    "step_index": step_index,
                    "template_id": result.template_id,
                    "operator_id": result.operator_id,
                    "available_precedent_count": result.available_precedent_count,
                    "returned_precedent_count": len(result.matches),
                    "match_evidence_ids": [
                        f"{prefix}.{index}"
                        for index in range(1, len(result.matches) + 1)
                    ],
                },
                provenance="deterministic_inference",
                schema_versions={"step_precedents": result.schema_version},
                uncertainty=(
                    "Library precedent supports an operator context; it does not "
                    "establish success for the planned substrate."
                ),
            )
        ]
        evidence.extend(
            EvidenceItem(
                evidence_id=f"{prefix}.{index}",
                layer="observation",
                source_id=match.match_id,
                payload_type="admitted_step_precedent",
                payload={
                    "route_alias": alias,
                    "step_index": step_index,
                    **match.to_dict(),
                },
                provenance="observed",
                schema_versions={"step_precedents": result.schema_version},
                uncertainty=(
                    "This is an admitted source observation, not an experimental "
                    "validation of the planned step."
                ),
            )
            for index, match in enumerate(result.matches, start=1)
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=tuple(evidence),
            packet={
                "result_ref": result_ref,
                "route_alias": alias,
                "step_index": step_index,
                "evidence": [_evidence_packet(item) for item in evidence],
            },
            authoritative_result=self._responses[result_ref],
            register_result_ref=False,
        )

    def verify_route(
        self,
        result_ref: str,
        alias: str,
    ) -> CapabilityResult:
        """Run independent whole-route verification over an existing route."""

        route = self._route(result_ref, alias)
        report = verify_planned_route(route)
        evidence = EvidenceItem(
            evidence_id=f"{alias}.verification",
            layer="route",
            source_id=route.route_id,
            payload_type="route_verification",
            payload={"route_alias": alias, **report.to_dict()},
            provenance="deterministic_inference",
            schema_versions={"route_verification": report.schema_version},
            uncertainty=(
                "Deterministic route verification is not experimental proof of "
                "yield or operational feasibility."
            ),
        )
        return CapabilityResult(
            result_ref=result_ref,
            evidence=(evidence,),
            packet={
                "result_ref": result_ref,
                "route_alias": alias,
                "evidence": [_evidence_packet(evidence)],
            },
            authoritative_result=self._responses[result_ref],
            register_result_ref=False,
        )

    def inspect_route(
        self,
        result_ref: str,
        alias: str,
        *,
        step_index: int | None = None,
    ) -> CapabilityResult:
        projection = self._projection(result_ref)
        packet = projection.inspection_packet(alias, step_index)
        evidence_ids = {item["evidence_id"] for item in packet["evidence"]}
        return CapabilityResult(
            result_ref=result_ref,
            evidence=tuple(
                item for item in projection.evidence if item.evidence_id in evidence_ids
            ),
            packet=packet,
            authoritative_result=self._responses[result_ref],
        )

    def compare_routes(
        self,
        result_ref: str,
        aliases: Tuple[str, ...],
    ) -> CapabilityResult:
        projection = self._projection(result_ref)
        packet = projection.comparison_packet(aliases)
        evidence_ids = {item["evidence_id"] for item in packet["evidence"]}
        return CapabilityResult(
            result_ref=result_ref,
            evidence=tuple(
                item for item in projection.evidence if item.evidence_id in evidence_ids
            ),
            packet=packet,
            authoritative_result=self._responses[result_ref],
        )

    def result(self, result_ref: str) -> MultistepRetrosynthesisResponse:
        try:
            return self._responses[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale multistep result: {result_ref}") from exc

    def _plan(
        self,
        request: MultistepRetrosynthesisRequest,
        *,
        candidate_exclusions: tuple[RouteCandidateExclusion, ...] = (),
        accepted_route_ids: frozenset[str] = frozenset(),
    ) -> CapabilityResult:
        response = self._coworker.plan(
            request,
            candidate_exclusions=candidate_exclusions,
        )
        projection = project_multistep_response(response)
        self._store(
            request,
            response,
            projection,
            candidate_exclusions=candidate_exclusions,
            accepted_route_ids=accepted_route_ids,
        )
        initial_ids = {
            item["evidence_id"] for item in projection.initial_packet()["evidence"]
        }
        return CapabilityResult(
            result_ref=projection.result_ref,
            evidence=tuple(
                item for item in projection.evidence if item.evidence_id in initial_ids
            ),
            packet=projection.initial_packet(),
            authoritative_result=response,
        )

    def _store(
        self,
        request: MultistepRetrosynthesisRequest,
        response: MultistepRetrosynthesisResponse,
        projection: MultistepEvidenceProjection,
        *,
        candidate_exclusions: tuple[RouteCandidateExclusion, ...] = (),
        accepted_route_ids: frozenset[str] = frozenset(),
    ) -> None:
        self._responses[projection.result_ref] = response
        self._requests[projection.result_ref] = request
        self._projections[projection.result_ref] = projection
        self._candidate_exclusions[projection.result_ref] = candidate_exclusions
        self._accepted_route_ids[projection.result_ref] = accepted_route_ids

    def _route(
        self,
        result_ref: str,
        alias: str,
    ) -> MultistepRetrosynthesisRoute:
        projection = self._projection(result_ref)
        try:
            route_id = dict(projection.route_aliases)[alias]
        except KeyError as exc:
            raise ValueError(f"unknown multistep route alias: {alias}") from exc
        response = self.result(result_ref)
        if response.result is None:
            raise ValueError("multistep response has no route result")
        routes = (*response.result.routes, *response.result.partial_routes)
        route = next((item for item in routes if item.route_id == route_id), None)
        if route is None:
            raise ValueError("multistep route alias resolved to a missing route")
        return route

    def _projection(self, result_ref: str) -> MultistepEvidenceProjection:
        try:
            return self._projections[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale multistep result: {result_ref}") from exc

    def _request(self, result_ref: str) -> MultistepRetrosynthesisRequest:
        try:
            return self._requests[result_ref]
        except KeyError as exc:
            raise ValueError(f"unknown or stale multistep result: {result_ref}") from exc


def _packet(
    result_ref: str,
    evidence: Tuple[EvidenceItem, ...],
    **aliases: str,
) -> Dict[str, Any]:
    return {
        "result_ref": result_ref,
        **aliases,
        "evidence": [
            {
                "evidence_id": item.evidence_id,
                "payload_type": item.payload_type,
                "payload": dict(item.payload),
                "uncertainty": item.uncertainty,
            }
            for item in evidence
        ],
    }


def _evidence_packet(item: EvidenceItem) -> Dict[str, Any]:
    """Return the compact evidence representation used by capability packets."""

    return {
        "evidence_id": item.evidence_id,
        "layer": item.layer,
        "payload_type": item.payload_type,
        "payload": dict(item.payload),
        "provenance": item.provenance,
        "uncertainty": item.uncertainty,
    }
