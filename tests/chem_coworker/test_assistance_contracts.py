"""Contract and transition tests for bounded assistance state."""

from __future__ import annotations

from dataclasses import replace

import pytest

from chem_coworker.assistance import (
    ActionRecord,
    AdvisoryClaim,
    AssistanceRequest,
    ConfirmedConstraint,
    EvidenceItem,
    load_assistance_policy,
    new_session,
    session_from_json,
)
from chem_coworker.assistance.contracts import (
    add_evidence,
    finish_session,
    record_action,
)


def _request() -> AssistanceRequest:
    return AssistanceRequest(
        objective="Explain and compare the recommended conditions",
        mode="conditions",
        structure_input="CCBr.O>>CCO",
    )


def _evidence() -> EvidenceItem:
    return EvidenceItem(
        evidence_id="candidate-1.ranking",
        layer="recommendation",
        source_id="result-1",
        payload_type="recommendation_score_trace",
        payload={"rank": 1, "score": 0.9},
        provenance="deterministic_inference",
    )


def test_session_ids_and_serialization_are_deterministic() -> None:
    first = new_session(_request())
    second = new_session(_request())

    assert first.session_id == second.session_id
    assert first.to_json() == second.to_json()
    assert first.to_dict()["status"] == "active"
    assert session_from_json(first.to_json()).to_json() == first.to_json()


def test_model_proposed_constraint_cannot_become_hard_filter() -> None:
    with pytest.raises(ValueError, match="cannot be hard"):
        ConfirmedConstraint(
            constraint_id="constraint-1",
            owner="condition_registry",
            kind="excluded_substance",
            value="substance:water",
            provenance="model_proposed",
            hard=True,
        )


def test_actions_are_evidence_checked_and_cannot_repeat() -> None:
    state = add_evidence(new_session(_request()), (_evidence(),))
    action = ActionRecord(
        action_id="action-1",
        action_name="inspect_condition_candidate",
        normalized_arguments={"candidate_alias": "candidate-1"},
        cited_evidence_ids=("candidate-1.ranking",),
        status="completed",
        provider_attempts=1,
    )
    advanced = record_action(state, action)

    assert advanced.turn == 1
    assert advanced.usage.provider_attempts == 1
    with pytest.raises(ValueError, match="identical action"):
        record_action(advanced, replace(action, action_id="action-2"))


def test_terminal_claims_must_reference_known_evidence() -> None:
    state = add_evidence(new_session(_request()), (_evidence(),))
    claim = AdvisoryClaim(
        claim_type="ranking_summary",
        subject_id="candidate-1",
        severity="info",
        statement="Candidate 1 is deterministically ranked first.",
        evidence_ids=("candidate-1.ranking",),
        uncertainty="The ranking is conditional on indexed precedents.",
    )
    finished = finish_session(
        state,
        status="completed",
        stopping_reason="Requested comparison is supported by available evidence.",
        claims=(claim,),
    )
    assert finished.is_terminal

    invalid = replace(claim, evidence_ids=("missing",))
    with pytest.raises(ValueError, match="unknown evidence"):
        finish_session(
            state,
            status="completed",
            stopping_reason="invalid",
            claims=(invalid,),
        )


def test_versioned_policy_is_closed_and_validates_bounds() -> None:
    policy = load_assistance_policy()

    assert policy.actions_for("conditions")[-1] == "finish"
    assert "retry_route_search" not in policy.actions_for("conditions")
    assert policy.default_rollout_state == "off"
    policy.validate_budget(policy.default_budget)
