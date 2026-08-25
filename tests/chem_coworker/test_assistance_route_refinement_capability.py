"""End-to-end capability test for deterministic route refinement."""

from __future__ import annotations

from dataclasses import replace

from cas_tools import molecule_identity
from chem_coworker.assistance import AssistanceRequest
from chem_coworker.assistance.capabilities import MultistepCapabilities
from chem_coworker.contracts import MultistepRetrosynthesisResponse
from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    RetrosynthesisConditionEvidence,
    plan_multistep_routes,
)


class _EmptyStock:
    def lookup(self, identity, *, provenance_limit: int = 5):
        return None


def _candidate(precursors: str, identity: str, score: float):
    candidate = GenericDisconnectionCandidate(
        target_smiles="CCCCCCCC",
        precursor_smiles=precursors,
        proposed_reaction_smiles=f"{precursors}>>CCCCCCCC",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id=f"template:{identity}",
        score=score,
        context_similarity=0.0,
        product_similarity=score,
        precursor_similarity=score,
        template_specificity=1.0,
        independent_reference_support=1,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=f"site:{identity}",
        precedent_reaction_ids=("reaction",),
        operator_id=f"operator:{identity}",
        realization_id=f"realization:{identity}",
        operator_signature=f"signature:{identity}",
        synthon_signature=f"synthon:{identity}",
    )
    assert molecule_identity(candidate.precursor_smiles) is not None
    return candidate


def _condition(reaction_smiles: str) -> RetrosynthesisConditionEvidence:
    insufficient = reaction_smiles.startswith("C.C>>")
    return RetrosynthesisConditionEvidence(
        status="insufficient_evidence" if insufficient else "recommended_direct",
        query_reaction_smiles=reaction_smiles,
        recommender_valid=not insufficient,
        recommendation_mode="verified_signature",
        retrieval_level=None if insufficient else "L2",
        uses_type_agnostic_fallback=False,
        candidate_count=0 if insufficient else 1,
        independent_candidate_count=0 if insufficient else 1,
        compatible_candidate_count=0 if insufficient else 1,
        independent_compatible_candidate_count=0 if insufficient else 1,
        excluded_candidate_count=0,
        best_recipe_score=None if insufficient else 0.8,
        best_recipe_compatibility_score=None if insufficient else 1.0,
        best_recipe_reference_support=0 if insufficient else 1,
        recommendations=(),
        warnings=("no condition precedent",) if insufficient else (),
        error=None,
    )


class _DeterministicCoworker:
    def __init__(self) -> None:
        self.exclusions = ()
        self.source = _candidate("C.C", "source", 0.99)
        self.alternative = _candidate("N.O", "alternative", 0.90)

    def plan(self, request, *, candidate_exclusions=()):
        self.exclusions = candidate_exclusions

        def expand(product: str, top_k: int):
            if product != "CCCCCCCC":
                return ()
            return (self.source, self.alternative)[:top_k]

        result = plan_multistep_routes(
            request.target_smiles,
            object(),
            _EmptyStock(),
            max_depth=request.max_depth,
            molecular_weight_threshold=50.0,
            top_k_routes=request.top_k,
            beam_width=request.beam_width,
            max_expansions=request.max_expansions,
            condition_evidence_evaluator=_condition,
            candidate_exclusions=candidate_exclusions,
            expander=expand,
        )
        return MultistepRetrosynthesisResponse(
            request=replace(request, molecular_weight_threshold=50.0),
            valid=True,
            result=result,
        )


def test_refinement_capability_resolves_aliases_and_preserves_lineage() -> None:
    coworker = _DeterministicCoworker()
    capabilities = MultistepCapabilities(coworker)  # type: ignore[arg-type]
    request = AssistanceRequest(
        objective="Resolve a condition gap using deterministic alternatives",
        mode="multistep",
        structure_input="CCCCCCCC",
    )
    initial = capabilities.plan_routes(request)
    response = capabilities.result(initial.result_ref)
    source_route = next(
        route
        for route in response.result.routes  # type: ignore[union-attr]
        if route.steps[0].candidate.strategy_id == coworker.source.strategy_id
    )
    projection = capabilities._projection(initial.result_ref)
    source_alias = next(
        alias
        for alias, route_id in projection.route_aliases
        if route_id == source_route.route_id
    )
    inspection = capabilities.inspect_route(
        initial.result_ref,
        source_alias,
        step_index=1,
    )
    issue = next(
        item
        for item in inspection.evidence
        if item.payload_type == "route_refinement_issue"
    )

    refined = capabilities.refine_route(
        request,
        initial.result_ref,
        route_alias=source_alias,
        step_index=1,
        objective="resolve_condition_gap",
        method="alternate_disconnection",
        issue_evidence_ids=(issue.evidence_id,),
        maximum_added_steps=0,
    )

    outcome = next(
        item
        for item in refined.evidence
        if item.payload_type == "route_refinement_outcome"
    )
    assert outcome.payload["status"] == "improved_alternative_found"
    assert outcome.payload["source_result_ref"] == initial.result_ref
    assert outcome.payload["source_route_preserved"] is True
    assert coworker.exclusions[0].strategy_id == coworker.source.strategy_id
    assert capabilities.result(initial.result_ref) is response
