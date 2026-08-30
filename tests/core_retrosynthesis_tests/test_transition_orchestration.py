"""Deterministic transition-provider orchestration regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

from core_retrosynthesis import (
    CallableTransitionProvider,
    ExpandLeafAction,
    ExpansionLeaf,
    ExpansionState,
    GenericDisconnectionCandidate,
    TransitionProviderBatch,
    TransitionProviderMetadata,
    TransitionProviderOrchestrator,
)


def _candidate(
    target: str = "CCN",
    precursors: str = "CCBr.N",
    *,
    status: str = "verified_signature",
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles=target,
        precursor_smiles=precursors,
        proposed_reaction_smiles=f"{precursors}>>{target}",
        transformation_kind=None,
        abstraction_level="L2",
        compiler_engine="reaction_core",
        template_id=f"template:{precursors}",
        score=0.9,
        context_similarity=0.0,
        product_similarity=0.9,
        precursor_similarity=0.9,
        template_specificity=1.0,
        independent_reference_support=1,
        forward_validation_status=status,
        center_transition_key="center",
        disconnection_site_key=f"site:{precursors}",
        precedent_reaction_ids=("reaction",),
        operator_id="operator",
        realization_id=f"realization:{precursors}",
        operator_signature="signature",
        synthon_signature="synthon",
    )


def _provider(
    provider_id: str,
    candidates: tuple[GenericDisconnectionCandidate, ...],
    *,
    maximum_budget: int = 5,
) -> CallableTransitionProvider:
    return CallableTransitionProvider(
        metadata=TransitionProviderMetadata(
            provider_id=provider_id,
            display_name=provider_id,
            capability_tags=("test",),
            maximum_proposal_budget=maximum_budget,
        ),
        expansion_function=lambda _target, budget: TransitionProviderBatch(
            candidates=candidates[:budget],
            diagnostics={"generated": len(candidates)},
        ),
    )


def _state(target: str = "CCN") -> ExpansionState:
    return ExpansionState(
        state_id="state:1",
        leaves=(ExpansionLeaf("leaf:1", target),),
    )


def test_action_contains_only_stable_ids_and_budget() -> None:
    action = ExpandLeafAction(
        leaf_id="leaf:1",
        provider_id="provider:a",
        proposal_budget=3,
    )

    assert action.action_id
    assert set(action.__dict__) == {
        "leaf_id",
        "provider_id",
        "proposal_budget",
        "action_id",
    }
    assert "CCN" not in action.__dict__.values()
    assert action == ExpandLeafAction(
        leaf_id="leaf:1",
        provider_id="provider:a",
        proposal_budget=3,
    )


def test_orchestrator_admits_verified_candidates_and_deduplicates() -> None:
    candidate = _candidate()
    orchestrator = TransitionProviderOrchestrator(
        (_provider("provider:a", (candidate, candidate)),)
    )

    outcome = orchestrator.expand(
        _state(),
        ExpandLeafAction("leaf:1", "provider:a", 2),
    )

    assert outcome.raw_candidate_count == 2
    assert len(outcome.admitted) == 1
    assert outcome.admitted[0].provider_id == "provider:a"
    assert outcome.admitted[0].provider_rank == 1
    assert outcome.duplicate_candidate_count == 1
    assert outcome.rejected == ()
    assert outcome.provider_diagnostics == {"generated": 2}
    assert outcome.to_dict()["admitted"][0]["transition_id"]


def test_orchestrator_rechecks_provider_claims_and_compatibility_rejections() -> None:
    not_verified = _candidate(status="core_only")
    incompatible = replace(
        _candidate(precursors="CCCl.N"),
        reaction_compatibility_disposition="reject",
    )
    false_signature_claim = _candidate(target="COC", precursors="CBr.CO")
    provider = _provider(
        "provider:a",
        (not_verified, incompatible, false_signature_claim),
    )
    orchestrator = TransitionProviderOrchestrator((provider,))

    outcome = orchestrator.expand(
        _state("COC"),
        ExpandLeafAction("leaf:1", "provider:a", 3),
    )

    assert outcome.admitted == ()
    assert tuple(item.reason for item in outcome.rejected) == (
        "provider_candidate_not_signature_verified",
        "provider_candidate_target_mismatch",
        "provider_candidate_signature_not_reproduced",
    )


def test_mapped_validation_reaction_resolves_ambiguous_display_reaction() -> None:
    candidate = replace(
        _candidate(target="COC", precursors="CBr.CO"),
        condition_query_reaction_smiles=(
            "[CH3:1][Br:2].[CH3:3][OH:4]>>[CH3:3][O:4][CH3:1]"
        ),
    )
    orchestrator = TransitionProviderOrchestrator(
        (_provider("provider:a", (candidate,)),)
    )

    outcome = orchestrator.expand(
        _state("COC"),
        ExpandLeafAction("leaf:1", "provider:a", 1),
    )

    assert len(outcome.admitted) == 1
    assert outcome.rejected == ()


def test_mapped_validation_reaction_must_match_display_structures() -> None:
    candidate = replace(
        _candidate(target="COC", precursors="CBr.CO"),
        condition_query_reaction_smiles=(
            "[CH3:1][Br:2].[NH3:3]>>[CH3:1][NH2:3]"
        ),
    )
    orchestrator = TransitionProviderOrchestrator(
        (_provider("provider:a", (candidate,)),)
    )

    outcome = orchestrator.expand(
        _state("COC"),
        ExpandLeafAction("leaf:1", "provider:a", 1),
    )

    assert outcome.admitted == ()
    assert tuple(item.reason for item in outcome.rejected) == (
        "provider_candidate_validation_reaction_mismatch",
    )


def test_shadow_comparison_preserves_overlap_and_provider_ranks() -> None:
    shared = _candidate()
    alternative = _candidate(precursors="CCCl.N")
    orchestrator = TransitionProviderOrchestrator(
        (
            _provider("provider:a", (shared,)),
            _provider("provider:b", (alternative, shared)),
        )
    )

    report = orchestrator.compare_shadow(
        _state(),
        leaf_id="leaf:1",
        provider_ids=("provider:a", "provider:b"),
        proposal_budget=2,
    )

    assert len(report.outcomes) == 2
    assert len(report.unique_transitions) == 2
    assert report.shared_transition_count == 1
    shared_result = next(
        item
        for item in report.unique_transitions
        if item.candidate.precursor_smiles == "CCBr.N"
    )
    assert shared_result.provider_ranks == (
        ("provider:a", 1),
        ("provider:b", 2),
    )
    assert report.to_dict()["shared_transition_count"] == 1


def test_orchestrator_rejects_unknown_or_out_of_budget_actions() -> None:
    orchestrator = TransitionProviderOrchestrator(
        (_provider("provider:a", (_candidate(),), maximum_budget=1),)
    )

    with pytest.raises(ValueError, match="unknown expansion leaf"):
        orchestrator.expand(
            _state(),
            ExpandLeafAction("leaf:missing", "provider:a", 1),
        )
    with pytest.raises(ValueError, match="unknown transition provider"):
        orchestrator.expand(
            _state(),
            ExpandLeafAction("leaf:1", "provider:missing", 1),
        )
    with pytest.raises(ValueError, match="exceeds provider capability"):
        orchestrator.expand(
            _state(),
            ExpandLeafAction("leaf:1", "provider:a", 2),
        )


def test_nonexpandable_leaf_cannot_be_sent_to_a_provider() -> None:
    state = ExpansionState(
        state_id="state:1",
        leaves=(ExpansionLeaf("leaf:1", "CCN", expandable=False),),
    )
    orchestrator = TransitionProviderOrchestrator(
        (_provider("provider:a", (_candidate(),)),)
    )

    with pytest.raises(ValueError, match="not expandable"):
        orchestrator.expand(
            state,
            ExpandLeafAction("leaf:1", "provider:a", 1),
        )
