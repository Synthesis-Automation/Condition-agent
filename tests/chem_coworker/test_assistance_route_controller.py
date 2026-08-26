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
        self.refinement = None

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

    def inspect_route(self, result_ref, alias, *, step_index=None):
        assert result_ref == "route-result-1"
        assert alias == "route-1"
        assert step_index == 1
        return CapabilityResult(
            result_ref=result_ref,
            evidence=(
                _evidence("route-1.step-1", {"step_index": 1}),
                EvidenceItem(
                    evidence_id="route-1.issue-1",
                    layer="route",
                    source_id=result_ref,
                    payload_type="route_refinement_issue",
                    payload={
                        "route_alias": "route-1",
                        "issue_id": "RISS1:test",
                        "kind": "condition_gap",
                        "step_index": 1,
                    },
                    provenance="deterministic_inference",
                ),
                EvidenceItem(
                    evidence_id="route-1.issue-1.repair-2",
                    layer="route",
                    source_id=result_ref,
                    payload_type="route_repair_proposal",
                    payload={
                        "route_alias": "route-1",
                        "step_index": 1,
                        "proposal_id": "RPROP1:unavailable",
                        "issue_id": "RISS1:test",
                        "status": "unavailable",
                        "objective": "resolve_condition_gap",
                        "refinement_method": None,
                        "maximum_added_steps": 2,
                    },
                    provenance="deterministic_inference",
                ),
                EvidenceItem(
                    evidence_id="route-1.issue-1.repair-1",
                    layer="route",
                    source_id=result_ref,
                    payload_type="route_repair_proposal",
                    payload={
                        "route_alias": "route-1",
                        "step_index": 1,
                        "proposal_id": "RPROP1:actionable",
                        "issue_id": "RISS1:test",
                        "status": "actionable",
                        "objective": "resolve_condition_gap",
                        "refinement_method": "alternate_disconnection",
                        "maximum_added_steps": 1,
                    },
                    provenance="deterministic_inference",
                ),
                EvidenceItem(
                    evidence_id="route-1.issue-1.repair-3",
                    layer="route",
                    source_id=result_ref,
                    payload_type="route_repair_proposal",
                    payload={
                        "route_alias": "route-1",
                        "step_index": 1,
                        "proposal_id": "RPROP1:second-actionable",
                        "issue_id": "RISS1:test",
                        "status": "actionable",
                        "objective": "resolve_condition_gap",
                        "refinement_method": "alternate_realization",
                        "maximum_added_steps": 1,
                    },
                    provenance="deterministic_inference",
                ),
            ),
            packet={},
            authoritative_result=self.first,
        )

    def search_step_precedents(self, result_ref, alias, *, step_index):
        assert result_ref == "route-result-1"
        assert alias == "route-1"
        assert step_index == 1
        return CapabilityResult(
            result_ref=result_ref,
            evidence=(
                _evidence(
                    "route-1.step-1.precedents.summary",
                    {"returned_precedent_count": 1},
                ),
            ),
            packet={},
            authoritative_result=self.first,
            register_result_ref=False,
        )

    def verify_route(self, result_ref, alias):
        assert result_ref == "route-result-1"
        assert alias == "route-1"
        return CapabilityResult(
            result_ref=result_ref,
            evidence=(
                _evidence(
                    "route-1.verification",
                    {"status": "verified_with_cautions"},
                ),
            ),
            packet={},
            authoritative_result=self.first,
            register_result_ref=False,
        )

    def apply_repair(self, request, result_ref, **kwargs):
        assert result_ref == "route-result-1"
        self.refinement = kwargs
        return CapabilityResult(
            result_ref="route-result-2",
            evidence=(
                _evidence(
                    "refinement.RINT1:test",
                    {
                        "status": "improved_alternative_found",
                        "source_route_preserved": True,
                    },
                ),
                _evidence("route-2.summary", {"solved": True}),
            ),
            packet={},
            authoritative_result=self.second,
        )

    def result(self, result_ref):
        return self.second if result_ref == "route-result-2" else self.first


class _RollbackThenImproveCapabilities(_MultistepCapabilities):
    def __init__(self):
        super().__init__()
        self.attempts = []

    def apply_repair(self, request, result_ref, **kwargs):
        self.attempts.append(kwargs["proposal_id"])
        if kwargs["proposal_id"] == "RPROP1:actionable":
            return CapabilityResult(
                result_ref="route-result-1",
                evidence=(
                    _evidence(
                        "refinement.RINT1:rejected",
                        {
                            "status": "alternatives_found_no_verified_improvement",
                            "retained_result_ref": "route-result-1",
                        },
                    ),
                ),
                packet={},
                authoritative_result=self.first,
                register_result_ref=False,
            )
        result = super().apply_repair(request, result_ref, **kwargs)
        return CapabilityResult(
            result_ref=result.result_ref,
            evidence=(
                _evidence(
                    "refinement.RINT1:accepted",
                    {
                        "status": "improved_alternative_found",
                        "retained_result_ref": "route-result-2",
                    },
                ),
                _evidence("route-2.summary", {"solved": True}),
            ),
            packet=result.packet,
            authoritative_result=result.authoritative_result,
        )


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


def test_target_audit_adds_read_only_evidence_before_one_step_search() -> None:
    transport = _QueueTransport(
        _action("audit_target"),
        _action("disconnect_target", evidence=("target.audit",)),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "Target and strategy evidence are available.",
            },
            evidence=("target.audit", "strategy-1.summary"),
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
            objective="Audit then disconnect the target",
            mode="retro",
            structure_input="C[C@H](O)Cl",
        )
    )

    audit = next(
        item for item in run.state.evidence if item.payload_type == "target_audit"
    )
    assert run.state.status == "completed"
    assert audit.payload["valid"] is True
    assert audit.payload["stereocenters"][0]["assigned"] is True
    assert run.state.domain_result_refs == ("retro-result-1",)


def test_multistep_read_only_tools_do_not_replace_authoritative_route_result() -> None:
    transport = _QueueTransport(
        _action("plan_routes"),
        _action(
            "search_step_precedents",
            arguments={"route_alias": "route-1", "step_index": 1},
            evidence=("route-1.summary",),
        ),
        _action(
            "verify_route",
            arguments={"route_alias": "route-1"},
            evidence=("route-1.step-1.precedents.summary",),
        ),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "The deterministic verification is available.",
            },
            evidence=("route-1.verification",),
            claims=(
                _claim(
                    "route-1",
                    "route-1.verification",
                    "The route passed deterministic verification with cautions.",
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
            objective="Inspect precedents and verify the route",
            mode="multistep",
            structure_input="CCN",
        )
    )

    assert run.state.status == "completed"
    assert run.state.domain_result_refs == ("route-result-1",)
    assert run.authoritative_result is capabilities.first


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


def test_multistep_repair_uses_only_deterministic_proposal_identity() -> None:
    transport = _QueueTransport(
        _action("plan_routes"),
        _action(
            "inspect_route_step",
            arguments={"route_alias": "route-1", "step_index": 1},
            evidence=("route-1.summary",),
        ),
        _action(
            "apply_repair",
            arguments={"proposal_id": "RPROP1:actionable"},
            evidence=("route-1.issue-1.repair-1",),
        ),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "A deterministic alternative was generated.",
            },
            evidence=("refinement.RINT1:test",),
            claims=(
                _claim(
                    "route-refinement",
                    "refinement.RINT1:test",
                    "The deterministic refinement found an alternative.",
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
            objective="Resolve the route condition gap without editing chemistry",
            mode="multistep",
            structure_input="CCN",
        )
    )

    assert run.state.status == "completed"
    assert run.state.usage.search_expansions == 1
    assert capabilities.refinement == {"proposal_id": "RPROP1:actionable"}
    refine_action = run.state.action_history[2]
    assert "smiles" not in str(refine_action.normalized_arguments).casefold()


def test_multistep_repair_can_retry_after_deterministic_rollback() -> None:
    transport = _QueueTransport(
        _action("plan_routes"),
        _action(
            "inspect_route_step",
            arguments={"route_alias": "route-1", "step_index": 1},
            evidence=("route-1.summary",),
        ),
        _action(
            "apply_repair",
            arguments={"proposal_id": "RPROP1:actionable"},
            evidence=("route-1.issue-1.repair-1",),
        ),
        _action(
            "apply_repair",
            arguments={"proposal_id": "RPROP1:second-actionable"},
            evidence=("route-1.issue-1.repair-3",),
        ),
        _action(
            "finish",
            arguments={
                "terminal_status": "completed",
                "stopping_reason": "The second deterministic repair improved the route.",
            },
            evidence=("refinement.RINT1:accepted",),
            claims=(
                _claim(
                    "route-refinement",
                    "refinement.RINT1:accepted",
                    "The second bounded repair was accepted.",
                ),
            ),
        ),
    )
    capabilities = _RollbackThenImproveCapabilities()
    controller = AssistanceController(
        transport=transport,
        multistep_capabilities=capabilities,  # type: ignore[arg-type]
    )

    run = controller.run(
        AssistanceRequest(
            objective="Try another deterministic proposal after rollback",
            mode="multistep",
            structure_input="CCN",
        )
    )

    assert run.state.status == "completed"
    assert run.state.usage.search_expansions == 2
    assert run.state.domain_result_refs == ("route-result-1", "route-result-2")
    assert capabilities.attempts == [
        "RPROP1:actionable",
        "RPROP1:second-actionable",
    ]


def test_multistep_repair_rejects_unavailable_proposal() -> None:
    transport = _QueueTransport(
        _action("plan_routes"),
        _action(
            "inspect_route_step",
            arguments={"route_alias": "route-1", "step_index": 1},
            evidence=("route-1.summary",),
        ),
        _action(
            "apply_repair",
            arguments={"proposal_id": "RPROP1:unavailable"},
            evidence=("route-1.issue-1.repair-2",),
        ),
    )
    capabilities = _MultistepCapabilities()
    controller = AssistanceController(
        transport=transport,
        multistep_capabilities=capabilities,  # type: ignore[arg-type]
    )

    run = controller.run(
        AssistanceRequest(
            objective="Do not execute unavailable chemistry",
            mode="multistep",
            structure_input="CCN",
        )
    )

    assert run.state.status == "blocked_by_policy"
    assert "actionable deterministic repair proposal" in (
        run.state.stopping_reason or ""
    )
    assert capabilities.refinement is None
