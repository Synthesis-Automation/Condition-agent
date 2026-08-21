"""Shared controller dispatch tests for retrosynthesis assistance modes."""

from __future__ import annotations

from chem_coworker.assistance import (
    AssistanceController,
    AssistanceRequest,
    AssistanceTransportResult,
)
from chem_coworker.assistance.capabilities import CapabilityResult
from chem_coworker.assistance.contracts import EvidenceItem
from chem_coworker.assistance.transport import AssistanceActionPayload


def _evidence(evidence_id: str, payload: dict) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=evidence_id,
        layer="route",
        source_id="fake-result",
        payload_type="fixture",
        payload=payload,
        provenance="deterministic_inference",
    )


def _action(name, *, arguments=None, evidence=(), claims=()):
    return {
        "action_name": name,
        "arguments": arguments or {},
        "cited_evidence_ids": list(evidence),
        "decision_summary": f"Use {name}.",
        "claims": list(claims),
    }


def _claim(subject: str, evidence_id: str, statement: str) -> dict:
    return {
        "claim_type": "route_summary",
        "subject_id": subject,
        "severity": "info",
        "statement": statement,
        "evidence_ids": [evidence_id],
        "uncertainty": "This is advisory over deterministic evidence.",
        "suggested_user_action": None,
    }


class _QueueTransport:
    def __init__(self, *payloads):
        self.payloads = list(payloads)
        self.packets = []

    def complete(self, packet, settings, *, max_retries):
        self.packets.append(packet)
        return AssistanceTransportResult(
            payload=AssistanceActionPayload.model_validate(self.payloads.pop(0))
        )


class _RetroCapabilities:
    def __init__(self):
        self.value = object()

    def disconnect_target(self, request):
        return CapabilityResult(
            result_ref="retro-result-1",
            evidence=(
                _evidence(
                    "strategy-1.summary",
                    {"forward_validation_status": "verified_signature"},
                ),
            ),
            packet={},
            authoritative_result=self.value,
        )

    def result(self, result_ref):
        assert result_ref == "retro-result-1"
        return self.value


class _MultistepCapabilities:
    def __init__(self):
        self.first = object()
        self.second = object()
        self.delta = None

    def plan_routes(self, request):
        return CapabilityResult(
            result_ref="route-result-1",
            evidence=(
                _evidence("query.multistep", {"route_kind": "partial"}),
                _evidence("route-1.summary", {"solved": False}),
            ),
            packet={},
            authoritative_result=self.first,
        )

    def retry_route_search(self, request, result_ref, delta):
        assert result_ref == "route-result-1"
        self.delta = delta
        return CapabilityResult(
            result_ref="route-result-2",
            evidence=(
                _evidence("query.multistep.expanded", {"route_kind": "solved"}),
                _evidence("route-2.summary", {"solved": True}),
            ),
            packet={},
            authoritative_result=self.second,
        )

    def result(self, result_ref):
        return self.second if result_ref == "route-result-2" else self.first


def test_same_controller_completes_one_step_mode_without_structure_edits() -> None:
    transport = _QueueTransport(
        _action("disconnect_target"),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "Validated strategy is available.",
            },
            evidence=("strategy-1.summary",),
            claims=(
                _claim(
                    "strategy-1",
                    "strategy-1.summary",
                    "Strategy 1 is forward validated.",
                ),
            ),
        ),
    )
    capabilities = _RetroCapabilities()
    controller = AssistanceController(
        transport=transport,
        retrosynthesis_capabilities=capabilities,  # type: ignore[arg-type]
    )

    run = controller.run(
        AssistanceRequest(
            objective="Explain the best disconnection",
            mode="retro",
            structure_input="CCN",
        )
    )

    assert run.state.status == "completed"
    assert run.authoritative_result is capabilities.value
    assert [item.action_name for item in run.state.action_history] == [
        "disconnect_target",
        "finish",
    ]


def test_multistep_retry_is_one_bounded_policy_delta_and_preserves_old_result() -> None:
    transport = _QueueTransport(
        _action("plan_routes"),
        _action(
            "retry_route_search",
            arguments={"search_depth_delta": 1, "max_expansions_delta": 2},
            evidence=("query.multistep",),
        ),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "Expanded search returned a solved route.",
            },
            evidence=("route-2.summary",),
            claims=(
                _claim(
                    "route-2",
                    "route-2.summary",
                    "The expanded deterministic search returned a solved route.",
                ),
            ),
        ),
    )
    capabilities = _MultistepCapabilities()
    controller = AssistanceController(
        transport=transport,
        multistep_capabilities=capabilities,  # type: ignore[arg-type]
    )

    run = controller.run(
        AssistanceRequest(
            objective="Find and explain a bounded route",
            mode="multistep",
            structure_input="CCN",
        )
    )

    assert run.state.status == "completed"
    assert run.state.domain_result_refs == ("route-result-1", "route-result-2")
    assert run.state.usage.search_expansions == 1
    assert capabilities.delta.max_depth_delta == 1
    assert capabilities.delta.max_expansions_delta == 2
    assert run.authoritative_result is capabilities.second
