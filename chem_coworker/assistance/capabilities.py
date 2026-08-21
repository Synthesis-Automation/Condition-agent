"""Closed application capabilities over authoritative domain services."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Tuple

from condition_recommender.models import GenericRecommendationResult
from core_retrosynthesis import RouteSearchPolicyDelta, apply_route_search_delta

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
        self._responses[projection.result_ref] = response
        self._requests[projection.result_ref] = request
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
