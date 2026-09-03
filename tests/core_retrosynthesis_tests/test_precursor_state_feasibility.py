"""Precursor-state feasibility and route-readiness regressions."""

from __future__ import annotations

from dataclasses import replace

import pytest

from core_retrosynthesis import (
    LatentStateTransition,
    MoleculeOccurrenceNode,
    PrecursorStateFeasibility,
    StepPrecedentLookupResult,
    StepPrecedentMatch,
    aggregate_precursor_state_route_feasibility,
    assess_precursor_state_feasibility,
)


def _transition() -> LatentStateTransition:
    return LatentStateTransition(
        transition_id="transition:test",
        step_id="step:test",
        retrosynthetic_depth=1,
        input_state_ids=("state:input",),
        output_state_id="state:output",
        tactical_input_occurrence_ids=(),
        transition_kind="state_change",
        action_class="unary_state_change",
        operator_id="operator:test",
        evidence_status="validated_route_step_projection",
    )


def _exact_lookup() -> StepPrecedentLookupResult:
    return StepPrecedentLookupResult(
        step_id="step:test",
        template_id="template:test",
        operator_id="operator:test",
        matches=(
            StepPrecedentMatch(
                match_id="match:test",
                reaction_id="reaction:test",
                reference_id="reference:test",
                template_id="template:test",
                operator_id="operator:test",
                product_smiles="CCO",
                precursor_smiles="CCBr.O",
                mapped_reaction_smiles="[CH3:1][CH2:2][Br:3].[OH2:4]>>[CH3:1][CH2:2][OH:4]",
                product_similarity=1.0,
                precursor_similarity=1.0,
            ),
        ),
        available_precedent_count=1,
    )


def test_weak_terminal_input_holds_an_otherwise_supported_route() -> None:
    weak_leaf = MoleculeOccurrenceNode(
        occurrence_id="leaf:test",
        smiles="CCBr",
        depth=1,
        terminal=True,
        terminal_evidence="heuristic_molecular_weight",
        unresolved_reason=None,
    )
    step = assess_precursor_state_feasibility(
        step_id="step:test",
        transition=_transition(),
        precursor_nodes=(weak_leaf,),
        precedent_lookup=_exact_lookup(),
        structurally_valid=True,
        hard_incompatible=False,
        compatibility_dispositions=("pass",),
        condition_status="recommended_direct",
    )

    assert step.evidence_level == "E4"
    assert step.promotion_recommendation == "eligible_for_route_review"
    assert step.weak_terminal_input_count == 1
    route = aggregate_precursor_state_route_feasibility((step,))
    assert route.status == "supported_with_cautions"
    assert route.promotion_recommendation == "hold_for_evidence"
    assert "ROUTE_HAS_WEAK_TERMINAL_INPUT_EVIDENCE" in route.warnings


def test_missing_transition_holds_exact_precedent_evidence() -> None:
    observed_leaf = MoleculeOccurrenceNode(
        occurrence_id="leaf:observed",
        smiles="CCBr",
        depth=1,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    step = assess_precursor_state_feasibility(
        step_id="step:test",
        transition=None,
        precursor_nodes=(observed_leaf,),
        precedent_lookup=_exact_lookup(),
        structurally_valid=True,
        hard_incompatible=False,
        compatibility_dispositions=("pass",),
        condition_status="recommended_direct",
    )

    assert step.evidence_level == "E4"
    assert step.support_status == "insufficient_evidence"
    assert step.promotion_recommendation == "hold_for_evidence"
    assert "LATENT_STATE_TRANSITION_UNAVAILABLE" in step.warnings


def test_serialized_feasibility_rejects_invalid_evidence() -> None:
    step = assess_precursor_state_feasibility(
        step_id="step:test",
        transition=_transition(),
        precursor_nodes=(),
        precedent_lookup=None,
        structurally_valid=True,
        hard_incompatible=False,
        compatibility_dispositions=("pass",),
        condition_status="insufficient_evidence",
    )

    assert PrecursorStateFeasibility.from_dict(step.to_dict()) == step
    with pytest.raises(ValueError, match="similarities must be in"):
        replace(step, best_precursor_similarity=1.1)
