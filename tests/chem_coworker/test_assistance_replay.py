"""Replay evaluation and blind-review packet tests."""

from __future__ import annotations

from chem_coworker.assistance import (
    AdvisoryClaim,
    AssistanceEvaluationCase,
    AssistanceRequest,
    EvidenceItem,
    build_blind_review_packet,
    evaluate_assistance_run,
    new_session,
)
from chem_coworker.assistance.contracts import add_evidence, finish_session


def _completed_state():
    state = new_session(
        AssistanceRequest(
            objective="Explain the recommendation",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
        )
    )
    state = add_evidence(
        state,
        (
            EvidenceItem(
                evidence_id="candidate-1.ranking",
                layer="recommendation",
                source_id="result-1",
                payload_type="ranking",
                payload={"rank": 1},
                provenance="deterministic_inference",
            ),
        ),
    )
    return finish_session(
        state,
        status="completed",
        stopping_reason="Evidence supports the answer.",
        claims=(
            AdvisoryClaim(
                claim_type="summary",
                subject_id="candidate-1",
                severity="info",
                statement="Candidate 1 is ranked first.",
                evidence_ids=("candidate-1.ranking",),
                uncertainty="Coverage is limited to the current index.",
            ),
        ),
    )


def test_evaluation_measures_replay_evidence_and_mutation() -> None:
    state = _completed_state()
    result = evaluate_assistance_run(
        AssistanceEvaluationCase(
            case_id="condition-explanation-1",
            acceptable_statuses=("completed",),
            required_evidence_ids=("candidate-1.ranking",),
            forbidden_actions=("retry_route_search",),
        ),
        state,
        authoritative_before={"rank": 1},
        authoritative_after={"rank": 1},
    )

    assert result.passed
    assert result.required_evidence_coverage == 1.0
    assert result.deterministic_mutation_count == 0
    assert result.replay_equivalent


def test_blind_packet_is_hash_bound_and_has_blank_adjudication() -> None:
    first = build_blind_review_packet(
        _completed_state(),
        authoritative_result={"rank": 1},
        one_shot_baseline="Baseline answer",
    )
    second = build_blind_review_packet(
        _completed_state(),
        authoritative_result={"rank": 1},
        one_shot_baseline="Baseline answer",
    )

    assert first["packet_id"] == second["packet_id"]
    assert first["adjudication"]["correctness"] is None
    assert first["final_claims"][0]["evidence_ids"] == (
        "candidate-1.ranking",
    )
