"""One-step capability loop over existing deterministic strategies."""

from __future__ import annotations

from dataclasses import replace

from cas_tools import molecule_identity
from chem_coworker.assistance import AssistanceRequest
from chem_coworker.assistance.capabilities import RetrosynthesisCapabilities
from chem_coworker.contracts import (
    RetrosynthesisResponse,
    RetrosynthesisStrategyCondition,
)
from core_retrosynthesis import (
    GenericDisconnectionCandidate,
    RetrosynthesisConditionEvidence,
    detect_functional_group_competition,
    group_strategy_candidates,
)


SELECTIVITY_REACTION = "Cc1[nH]cnc1CCl.NCCS>>SCCNCc1c(C)[nH]cn1"
SELECTIVITY_PRODUCT = molecule_identity(
    SELECTIVITY_REACTION.split(">>", 1)[1]
).canonical_smiles
SELECTIVITY_PRECURSORS = SELECTIVITY_REACTION.split(">>", 1)[0]
SELECTIVITY_QUERY = f"{SELECTIVITY_PRECURSORS}>>{SELECTIVITY_PRODUCT}"


def _candidate(name: str, *, warning: bool) -> GenericDisconnectionCandidate:
    candidate = GenericDisconnectionCandidate(
        target_smiles="CCN",
        precursor_smiles="CCBr.N",
        proposed_reaction_smiles="CCBr.N>>CCN",
        transformation_kind="substitution",
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template:{name}",
        score=0.9,
        context_similarity=0.8,
        product_similarity=0.9,
        precursor_similarity=0.8,
        template_specificity=0.8,
        independent_reference_support=2,
        forward_validation_status="verified_signature",
        center_transition_key=f"center:{name}",
        disconnection_site_key=f"SITE1:{name}",
        precedent_reaction_ids=(f"precedent:{name}",),
        operator_id=f"OP1:{name}",
        realization_id=f"REAL1:{name}",
        operator_signature=f"operator:{name}",
        synthon_signature=f"SYN1:{name}",
    )
    if warning:
        competition = detect_functional_group_competition(SELECTIVITY_REACTION)
        assert competition is not None
        candidate = replace(candidate, selectivity_warnings=(competition,))
    return candidate


class _Coworker:
    def disconnect(self, request):
        source = group_strategy_candidates(
            (_candidate("source", warning=True),), top_k_strategies=1
        )[0]
        target = group_strategy_candidates(
            (_candidate("target", warning=False),), top_k_strategies=1
        )[0]
        return RetrosynthesisResponse(
            request=request,
            valid=True,
            strategies=(source, target),
        )

    def recommend_strategy_conditions(self, strategy, request):
        del strategy, request
        return None


def _selectivity_condition() -> RetrosynthesisConditionEvidence:
    recommendation = {
        "rank": 1,
        "recipe_id": "recipe:neutral",
        "recipe_core_id": "recipe-core:neutral",
        "resolved_recipe": {
            "medium": "neutral",
            "salt_state": "free_base",
            "temperature_c": 25.0,
        },
        "compatibility_score": 1.0,
        "reference_support": 2,
        "precedent_reaction_contexts": (
            {
                "reaction_smiles": SELECTIVITY_REACTION,
                "reference_id": "reference:1",
            },
            {
                "reaction_smiles": SELECTIVITY_REACTION,
                "reference_id": "reference:2",
            },
        ),
    }
    return RetrosynthesisConditionEvidence(
        status="recommended_direct",
        query_reaction_smiles=SELECTIVITY_QUERY,
        recommender_valid=True,
        recommendation_mode="verified_signature",
        retrieval_level="L2",
        uses_type_agnostic_fallback=False,
        candidate_count=2,
        independent_candidate_count=2,
        compatible_candidate_count=2,
        independent_compatible_candidate_count=2,
        excluded_candidate_count=0,
        best_recipe_score=0.9,
        best_recipe_compatibility_score=1.0,
        best_recipe_reference_support=2,
        recommendations=(recommendation,),
        warnings=(),
        error=None,
    )


class _SelectivityCoworker:
    def disconnect(self, request):
        warning = detect_functional_group_competition(SELECTIVITY_QUERY)
        assert warning is not None
        candidate = GenericDisconnectionCandidate(
            target_smiles=SELECTIVITY_PRODUCT,
            precursor_smiles=SELECTIVITY_PRECURSORS,
            proposed_reaction_smiles=SELECTIVITY_QUERY,
            transformation_kind="substitution",
            abstraction_level="L2",
            compiler_engine="test",
            template_id="template:selectivity",
            score=0.9,
            context_similarity=0.9,
            product_similarity=0.9,
            precursor_similarity=0.9,
            template_specificity=0.9,
            independent_reference_support=2,
            forward_validation_status="verified_signature",
            center_transition_key="center:selectivity",
            disconnection_site_key="SITE1:selectivity",
            precedent_reaction_ids=("reference:1", "reference:2"),
            operator_id="OP1:selectivity",
            realization_id="REAL1:selectivity",
            operator_signature="operator:selectivity",
            synthon_signature="SYN1:selectivity",
            selectivity_warnings=(warning,),
        )
        strategy = group_strategy_candidates(
            (candidate,), top_k_strategies=1
        )[0]
        return RetrosynthesisResponse(
            request=request,
            valid=True,
            strategies=(strategy,),
            condition_evidence=(
                RetrosynthesisStrategyCondition(
                    strategy_id=strategy.strategy_id,
                    evidence=_selectivity_condition(),
                ),
            ),
        )

    def recommend_strategy_conditions(self, strategy, request):
        del strategy, request
        raise AssertionError("condition repair must reuse the attached recipe")


def _request() -> AssistanceRequest:
    return AssistanceRequest(
        objective="Find and repair a one-step disconnection.",
        mode="retro",
        structure_input="CCN",
    )


def test_apply_repair_selects_existing_strategy_and_verifies_it() -> None:
    capabilities = RetrosynthesisCapabilities(_Coworker())  # type: ignore[arg-type]
    initial = capabilities.disconnect_target(_request())
    inspection = capabilities.inspect_strategy(initial.result_ref, "strategy-1")
    proposal = next(
        item
        for item in inspection.evidence
        if item.payload_type == "single_step_repair_proposal"
        and item.payload["status"] == "actionable"
        and item.payload["repair_kind"] == "alternate_strategy"
    )

    repaired = capabilities.apply_repair(
        _request(),
        initial.result_ref,
        proposal_id=str(proposal.payload["proposal_id"]),
    )

    assert repaired.result_ref != initial.result_ref
    assert repaired.authoritative_result.selected_strategy_id == (
        repaired.authoritative_result.strategies[1].strategy_id
    )
    outcome = next(
        item
        for item in repaired.evidence
        if item.payload_type == "single_step_refinement_outcome"
    )
    verification = next(
        item
        for item in repaired.evidence
        if item.payload_type == "single_step_refinement_verification"
    )
    assert outcome.payload["status"] == "improved_alternative_found"
    assert verification.payload["accepted"] is True
    assert verification.payload["status"] == "verified"


def test_verify_strategy_is_read_only() -> None:
    capabilities = RetrosynthesisCapabilities(_Coworker())  # type: ignore[arg-type]
    initial = capabilities.disconnect_target(_request())

    result = capabilities.verify_strategy(initial.result_ref, "strategy-2")

    assert result.result_ref == initial.result_ref
    assert result.register_result_ref is False
    assert result.evidence[0].payload["status"] == "verified"


def test_apply_condition_recipe_resolves_single_step_selectivity_issue() -> None:
    capabilities = RetrosynthesisCapabilities(
        _SelectivityCoworker()  # type: ignore[arg-type]
    )
    request = AssistanceRequest(
        objective="Resolve one-step selectivity from admitted condition evidence.",
        mode="retro",
        structure_input=SELECTIVITY_PRODUCT,
    )
    initial = capabilities.disconnect_target(request)
    inspection = capabilities.inspect_strategy(initial.result_ref, "strategy-1")
    proposal = next(
        item
        for item in inspection.evidence
        if item.payload_type == "single_step_repair_proposal"
        and item.payload["repair_kind"] == "condition_selectivity"
        and item.payload["status"] == "actionable"
    )

    repaired = capabilities.apply_repair(
        request,
        initial.result_ref,
        proposal_id=str(proposal.payload["proposal_id"]),
    )
    repaired_inspection = capabilities.inspect_strategy(
        repaired.result_ref, "strategy-1"
    )

    assert repaired.result_ref != initial.result_ref
    assert not any(
        item.payload_type == "single_step_refinement_issue"
        for item in repaired_inspection.evidence
    )
    condition = repaired.authoritative_result.condition_evidence[0]
    assert condition.condition_selectivity_assessment is not None
    assert condition.condition_selectivity_assessment.status == "supported"
