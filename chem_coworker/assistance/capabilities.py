"""Closed application capabilities over authoritative domain services."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Tuple

from condition_recommender.models import GenericRecommendationResult
from core_retrosynthesis import (
    MultistepRetrosynthesisRoute,
    RouteRefinementIntent,
    RouteSearchPolicyDelta,
    apply_route_search_delta,
    build_route_refinement_plan,
    collect_route_refinement_issues,
    summarize_route_refinement,
)

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
from chem_coworker.service import ConditionCoworker

from .contracts import AssistanceRequest, EvidenceItem, canonical_json
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
        return self._plan(expanded)

    def refine_route(
        self,
        request: AssistanceRequest,
        result_ref: str,
        *,
        route_alias: str,
        step_index: int,
        objective: str,
        method: str,
        issue_evidence_ids: Tuple[str, ...],
        maximum_added_steps: int,
    ) -> CapabilityResult:
        """Run one chemistry-directed re-search without model-authored structures."""

        source_route = self._route(result_ref, route_alias)
        if step_index < 1 or step_index > len(source_route.steps):
            raise ValueError("route refinement step index is out of range")
        projection = self._projection(result_ref)
        evidence_by_id = {item.evidence_id: item for item in projection.evidence}
        issue_ids = []
        for evidence_id in issue_evidence_ids:
            evidence = evidence_by_id.get(evidence_id)
            if evidence is None:
                raise ValueError(f"unknown route issue evidence: {evidence_id}")
            if evidence.payload_type != "route_refinement_issue":
                raise ValueError("route refinement requires typed issue evidence")
            if evidence.payload.get("route_alias") != route_alias:
                raise ValueError("route issue evidence belongs to another route")
            issue_id = str(evidence.payload.get("issue_id") or "")
            if not issue_id:
                raise ValueError("route issue evidence has no stable issue ID")
            issue_ids.append(issue_id)
        source_step = source_route.steps[step_index - 1]
        intent = RouteRefinementIntent(
            source_route_id=source_route.route_id,
            source_step_id=source_step.step_id,
            objective=objective,  # type: ignore[arg-type]
            method=method,  # type: ignore[arg-type]
            issue_ids=tuple(issue_ids),
            maximum_added_steps=maximum_added_steps,
        )
        source_issues = collect_route_refinement_issues(source_route)
        plan = build_route_refinement_plan(source_route, intent, source_issues)
        prior = self._request(result_ref)
        policy = request.route_search_policy
        depth = min(
            policy.maximum_max_depth,
            prior.max_depth + maximum_added_steps,
        )
        refined_request = MultistepRetrosynthesisRequest(
            **{
                **prior.__dict__,
                "max_depth": depth,
                "beam_width": policy.maximum_beam_width,
                "max_expansions": policy.maximum_max_expansions,
            }
        )
        response = self._coworker.plan(
            refined_request,
            candidate_exclusions=(plan.exclusion,),
        )
        refined_projection = project_multistep_response(response)
        if response.result is None:
            raise ValueError("deterministic refinement returned no route result")
        outcome = summarize_route_refinement(
            source_route,
            intent,
            response.result,
        )
        lineage = EvidenceItem(
            evidence_id=f"refinement.{intent.intent_id}",
            layer="route",
            source_id=refined_projection.result_ref,
            payload_type="route_refinement_outcome",
            payload={
                **outcome.to_dict(),
                "source_result_ref": result_ref,
                "refined_result_ref": refined_projection.result_ref,
                "excluded_candidate_scope": plan.excluded_candidate_scope,
            },
            provenance="deterministic_inference",
            schema_versions={"route_refinement": outcome.schema_version},
            uncertainty=(
                "A lower deterministic issue count is not experimental proof of "
                "reaction success."
            ),
        )
        self._store(refined_request, response, refined_projection)
        initial_ids = {
            item["evidence_id"]
            for item in refined_projection.initial_packet()["evidence"]
        }
        exposed = tuple(
            item
            for item in refined_projection.evidence
            if item.evidence_id in initial_ids
        )
        packet = {
            **refined_projection.initial_packet(),
            "refinement": lineage.payload,
        }
        return CapabilityResult(
            result_ref=refined_projection.result_ref,
            evidence=(lineage, *exposed),
            packet=packet,
            authoritative_result=response,
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

    def _plan(self, request: MultistepRetrosynthesisRequest) -> CapabilityResult:
        response = self._coworker.plan(request)
        projection = project_multistep_response(response)
        self._store(request, response, projection)
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
    ) -> None:
        self._responses[projection.result_ref] = response
        self._requests[projection.result_ref] = request
        self._projections[projection.result_ref] = projection

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
