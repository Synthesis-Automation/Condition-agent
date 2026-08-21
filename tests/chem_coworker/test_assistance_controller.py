"""Bounded controller tests using deterministic fake provider actions."""

from __future__ import annotations

from types import SimpleNamespace

from condition_recommender.models import (
    GenericConditionRecommendation,
    GenericRecommendationResult,
    RecommendationScoreTrace,
)

from chem_coworker.assistance import (
    AssistanceBudget,
    AssistanceController,
    AssistanceRequest,
    AssistanceTransportResult,
    ConditionCapabilities,
    EvidenceItem,
)
from chem_coworker.assistance.transport import AssistanceActionPayload


def _result() -> GenericRecommendationResult:
    score_trace = RecommendationScoreTrace(
        similarity_components={"edit": 1.0},
        similarity_contributions={"edit": 0.8},
        ranking_components={"similarity": 0.8},
        ranking_contributions={"similarity": 0.8},
        applied_ranking_weights={"similarity": 1.0},
        independent_evidence_count=2,
        observed_outcome_count=2,
        pool_yield_prior_pct=70.0,
        definition_versions={"ranking": "ranking.v1@1.0"},
    )
    candidate = GenericConditionRecommendation(
        rank=1,
        recipe_id="recipe:1",
        recipe_core_id="recipe-core:1",
        recipe_variant_ids=("variant:1",),
        resolved_recipe={"components": []},
        score=0.8,
        similarity_score=0.8,
        compatibility_score=1.0,
        expected_yield_pct=72.0,
        support=2,
        observation_support=2,
        reference_support=2,
        condition_series_support=1,
        dataset_support=1,
        retrieval_level="generic_signature",
        precedent_reaction_ids=("reaction:1",),
        precedent_reference_ids=("reference:1",),
        explanation=("Edit and environment match.",),
        score_trace=score_trace,
    )
    return GenericRecommendationResult(
        query_reaction_smiles="CCBr.O>>CCO",
        valid=True,
        query_signature_id="signature:1",
        retrieval_definition_version="retrieval.v1@1.0",
        retrieval_level="generic_signature",
        candidate_count=2,
        independent_candidate_count=2,
        compatible_candidate_count=1,
        independent_compatible_candidate_count=1,
        excluded_candidate_count=1,
        recommendations=(candidate,),
    )


class _FakeCoworker:
    def __init__(self) -> None:
        self.requests = []

    def recommend(self, request):
        self.requests.append(request)
        return SimpleNamespace(result=_result())


class _QueueTransport:
    def __init__(self, *payloads: dict) -> None:
        self.payloads = list(payloads)
        self.packets = []

    def complete(self, packet, settings, *, max_retries):
        self.packets.append(packet)
        return AssistanceTransportResult(
            payload=AssistanceActionPayload.model_validate(self.payloads.pop(0)),
            response_id=f"response-{len(self.packets)}",
            input_tokens=10,
            output_tokens=5,
            elapsed_ms=2,
        )


def _action(name: str, *, arguments=None, evidence=(), claims=()) -> dict:
    return {
        "action_name": name,
        "arguments": arguments or {},
        "cited_evidence_ids": list(evidence),
        "decision_summary": f"Use {name}.",
        "claims": list(claims),
    }


def test_controller_recommends_inspects_and_finishes_with_cited_claim() -> None:
    transport = _QueueTransport(
        _action("recommend_conditions"),
        _action(
            "inspect_condition_candidate",
            arguments={"candidate_alias": "candidate-1"},
            evidence=("candidate-1.summary",),
        ),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "The requested explanation is supported.",
            },
            evidence=("candidate-1.ranking",),
            claims=(
                {
                    "claim_type": "ranking_summary",
                    "subject_id": "candidate-1",
                    "severity": "info",
                    "statement": "Candidate 1 is the deterministic top result.",
                    "evidence_ids": ["candidate-1.ranking"],
                    "uncertainty": "The result depends on indexed precedent coverage.",
                    "suggested_user_action": None,
                },
            ),
        ),
    )
    coworker = _FakeCoworker()
    controller = AssistanceController(
        transport=transport,
        condition_capabilities=ConditionCapabilities(coworker),
    )

    run = controller.run(
        AssistanceRequest(
            objective="Explain the best condition recommendation",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
        )
    )

    assert run.state.status == "completed"
    assert run.authoritative_result is not None
    assert run.authoritative_result.to_dict() == _result().to_dict()
    assert [item.action_name for item in run.state.action_history] == [
        "recommend_conditions",
        "inspect_condition_candidate",
        "finish",
    ]
    assert run.state.usage.input_tokens == 30
    assert "candidate-1.ranking" in run.state.allowed_evidence_ids
    assert coworker.requests[0].review.mode == "off"


def test_controller_can_stop_for_one_material_clarification() -> None:
    transport = _QueueTransport(
        _action(
            "propose_clarification",
            arguments={
                "question_id": "question-1",
                "prompt": "Must palladium-containing recipes be excluded?",
                "constraint_owner": "condition_registry",
                "constraint_kind": "excluded_substance",
                "reason": "The answer would change the candidate set.",
            },
        )
    )
    controller = AssistanceController(
        transport=transport,
        condition_capabilities=ConditionCapabilities(_FakeCoworker()),
    )

    run = controller.run(
        AssistanceRequest(
            objective="Find conditions compatible with my constraints",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
        )
    )

    assert run.state.status == "needs_user_input"
    assert run.state.unresolved_questions[0].constraint_kind == "excluded_substance"
    assert run.state.usage.clarification_cycles == 1
    assert run.authoritative_result is None


def test_unknown_evidence_reference_fails_closed_before_capability_execution() -> None:
    coworker = _FakeCoworker()
    controller = AssistanceController(
        transport=_QueueTransport(
            _action(
                "recommend_conditions",
                evidence=("invented.evidence",),
            )
        ),
        condition_capabilities=ConditionCapabilities(coworker),
    )

    run = controller.run(
        AssistanceRequest(
            objective="Recommend conditions",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
        )
    )

    assert run.state.status == "blocked_by_policy"
    assert not coworker.requests


def test_provider_error_redacts_secret_like_values() -> None:
    class FailingTransport:
        def complete(self, packet, settings, *, max_retries):
            raise RuntimeError("credential sk-supersecret123456789 failed")

    controller = AssistanceController(
        transport=FailingTransport(),
        condition_capabilities=ConditionCapabilities(_FakeCoworker()),
    )
    run = controller.run(
        AssistanceRequest(
            objective="Recommend conditions",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
        )
    )

    assert run.state.status == "provider_failed"
    assert "supersecret" not in run.state.stopping_reason
    assert "[REDACTED]" in run.state.stopping_reason


def test_confirmed_clarification_is_normalized_then_recomputed_canonically() -> None:
    transport = _QueueTransport(
        _action(
            "propose_clarification",
            arguments={
                "question_id": "question-constraint",
                "prompt": "Exclude Pd(PPh3)4?",
                "constraint_owner": "condition_registry",
                "constraint_kind": "excluded_substance",
                "reason": "The answer can change the returned recipes.",
            },
        ),
        _action("recommend_conditions"),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "Recommendation was recomputed.",
            },
            evidence=("candidate-1.summary",),
            claims=(
                {
                    "claim_type": "constraint_application",
                    "subject_id": "candidate-1",
                    "severity": "info",
                    "statement": "The result was recomputed after confirmation.",
                    "evidence_ids": ["candidate-1.summary"],
                    "uncertainty": "Indexed precedent coverage remains limited.",
                    "suggested_user_action": None,
                },
            ),
        ),
    )
    coworker = _FakeCoworker()
    controller = AssistanceController(
        transport=transport,
        condition_capabilities=ConditionCapabilities(coworker),
    )
    request = AssistanceRequest(
        objective="Respect my catalyst restriction",
        mode="conditions",
        structure_input="CCBr.O>>CCO",
    )

    waiting = controller.run(request)
    run = controller.resume_with_condition_constraint(waiting.state, "Pd(PPh3)4")

    assert run.state.status == "completed"
    assert coworker.requests[0].condition_constraints.constraints[0].normalized_value == (
        "cas:14221-01-3"
    )
    assert any(
        item.evidence_id.startswith("user.constraint.CCON1:")
        for item in run.state.evidence
    )


def test_action_budget_exhaustion_is_a_typed_terminal_outcome() -> None:
    controller = AssistanceController(
        transport=_QueueTransport(_action("recommend_conditions")),
        condition_capabilities=ConditionCapabilities(_FakeCoworker()),
    )

    run = controller.run(
        AssistanceRequest(
            objective="Recommend conditions",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
            budget=AssistanceBudget(max_action_turns=1),
        )
    )

    assert run.state.status == "budget_exhausted"
    assert run.authoritative_result is not None
    assert run.state.stopping_reason == "maximum action turns exhausted"


def test_two_distinct_actions_without_new_evidence_stop_as_no_progress() -> None:
    summary = EvidenceItem(
        evidence_id="candidate-1.summary",
        layer="recommendation",
        source_id="result-1",
        payload_type="summary",
        payload={"rank": 1},
        provenance="deterministic_inference",
    )

    class NoProgressCapabilities:
        value = object()

        def request_fingerprint(self, request):
            return "request-1"

        def recommend_conditions(self, request):
            from chem_coworker.assistance.capabilities import CapabilityResult

            return CapabilityResult("result-1", (summary,), {}, self.value)

        def inspect_condition_candidate(self, **kwargs):
            from chem_coworker.assistance.capabilities import CapabilityResult

            return CapabilityResult("result-1", (summary,), {}, self.value)

        def compare_condition_candidates(self, **kwargs):
            from chem_coworker.assistance.capabilities import CapabilityResult

            return CapabilityResult("result-1", (summary,), {}, self.value)

        def result(self, result_ref):
            return self.value

    transport = _QueueTransport(
        _action("recommend_conditions"),
        _action(
            "inspect_condition_candidate",
            arguments={"candidate_alias": "candidate-1"},
            evidence=("candidate-1.summary",),
        ),
        _action(
            "compare_condition_candidates",
            arguments={"candidate_aliases": ["candidate-1", "candidate-2"]},
            evidence=("candidate-1.summary",),
        ),
    )
    controller = AssistanceController(
        transport=transport,
        condition_capabilities=NoProgressCapabilities(),  # type: ignore[arg-type]
    )

    run = controller.run(
        AssistanceRequest(
            objective="Compare conditions",
            mode="conditions",
            structure_input="CCBr.O>>CCO",
        )
    )

    assert run.state.status == "no_progress"
    assert len(run.state.action_history) == 3
