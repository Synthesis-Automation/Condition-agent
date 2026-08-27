"""Bounded state machine for evidence-grounded chem-coworker assistance."""

from __future__ import annotations

from dataclasses import dataclass, replace
import re
from typing import Any, Dict, Mapping, Tuple

from condition_registry import normalize_condition_constraint
from core_retrosynthesis import RouteSearchPolicyDelta

from .capabilities import (
    CapabilityResult,
    ChemistryCapabilities,
    ConditionCapabilities,
    MultistepCapabilities,
    RetrosynthesisCapabilities,
)
from .contracts import (
    ActionRecord,
    AdvisoryClaim,
    AssistanceRequest,
    AssistanceSessionState,
    ClarificationQuestion,
    EvidenceItem,
    add_evidence,
    finish_session,
    new_session,
    record_action,
    stable_assistance_id,
)
from .policy import AssistancePolicy, load_assistance_policy
from .transport import (
    AssistanceActionPayload,
    AssistanceProviderSettings,
    AssistanceTransport,
    AssistanceTransportResult,
)


_SECRET_PATTERN = re.compile(r"\b(?:sk|dashscope)-[A-Za-z0-9_-]{8,}\b")


@dataclass(frozen=True)
class AssistanceRunResult:
    """Terminal assistance state beside the last authoritative domain result."""

    state: AssistanceSessionState
    authoritative_result: object | None

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the advisory trace without altering the domain result."""

        return {
            "state": self.state.to_dict(),
            "authoritative_result": (
                self.authoritative_result.to_dict()
                if hasattr(self.authoritative_result, "to_dict")
                else self.authoritative_result
            ),
        }


class AssistanceController:
    """Choose, validate, execute, and record one closed action per turn."""

    def __init__(
        self,
        *,
        transport: AssistanceTransport,
        chemistry_capabilities: ChemistryCapabilities | None = None,
        condition_capabilities: ConditionCapabilities | None = None,
        retrosynthesis_capabilities: RetrosynthesisCapabilities | None = None,
        multistep_capabilities: MultistepCapabilities | None = None,
        policy: AssistancePolicy | None = None,
    ) -> None:
        self._transport = transport
        self._chemistry = chemistry_capabilities or ChemistryCapabilities()
        self._conditions = condition_capabilities
        self._retro = retrosynthesis_capabilities
        self._multistep = multistep_capabilities
        self._policy = policy or load_assistance_policy()

    def run(
        self,
        request: AssistanceRequest,
        *,
        initial_state: AssistanceSessionState | None = None,
    ) -> AssistanceRunResult:
        """Run until completion, clarification, failure, or a mandatory stop."""

        self._policy.validate_budget(request.budget)
        if self._capabilities_for(request.mode) is None:
            return AssistanceRunResult(
                state=finish_session(
                    initial_state or new_session(request),
                    status="blocked_by_policy",
                    stopping_reason=f"{request.mode} capabilities are not configured",
                ),
                authoritative_result=None,
            )
        state = initial_state or new_session(request)
        if state.request != request:
            raise ValueError("initial state request does not match controller request")
        if state.is_terminal:
            raise ValueError("cannot run a terminal assistance session")
        last_result = self._last_result(state)
        no_progress_turns = 0

        while not state.is_terminal:
            exhausted = self._budget_reason(state)
            if exhausted:
                state = finish_session(
                    state,
                    status="budget_exhausted",
                    stopping_reason=exhausted,
                )
                break
            allowed_actions = self._allowed_actions(state)
            packet = self._controller_packet(state, allowed_actions)
            settings = AssistanceProviderSettings.from_mapping(
                request.provider_settings
            )
            try:
                transport_result = self._transport.complete(
                    packet,
                    settings,
                    max_retries=request.budget.max_provider_retries_per_turn,
                )
            except Exception as exc:
                state = finish_session(
                    state,
                    status="provider_failed",
                    stopping_reason=self._safe_error("provider failure", exc),
                )
                break

            payload = transport_result.payload
            if (
                transport_result.attempts
                > 1 + request.budget.max_provider_retries_per_turn
            ):
                state = finish_session(
                    state,
                    status="blocked_by_policy",
                    stopping_reason="provider exceeded the per-turn retry allowance",
                )
                break
            if payload.action_name not in allowed_actions:
                state = finish_session(
                    state,
                    status="blocked_by_policy",
                    stopping_reason=(
                        f"action {payload.action_name!r} is not allowed in current state"
                    ),
                )
                break
            unknown = set(payload.cited_evidence_ids) - state.allowed_evidence_ids
            if unknown:
                state = finish_session(
                    state,
                    status="blocked_by_policy",
                    stopping_reason=f"action cited unknown evidence: {sorted(unknown)!r}",
                )
                break

            try:
                normalized_arguments = self._normalized_arguments(payload, state)
            except ValueError as exc:
                state = finish_session(
                    state,
                    status="blocked_by_policy",
                    stopping_reason=self._safe_error(
                        "invalid provider action",
                        exc,
                    ),
                )
                break
            action_id = stable_assistance_id(
                "ACTION",
                {
                    "session_id": state.session_id,
                    "turn": state.turn + 1,
                    "name": payload.action_name,
                    "arguments": normalized_arguments,
                },
            )
            projected_tokens = (
                state.usage.total_tokens
                + transport_result.input_tokens
                + transport_result.output_tokens
            )
            projected_elapsed = state.usage.elapsed_ms + transport_result.elapsed_ms
            if (
                projected_tokens > request.budget.max_total_tokens
                or projected_elapsed > request.budget.max_elapsed_ms
            ):
                action = self._action_record(
                    action_id,
                    payload,
                    normalized_arguments,
                    transport_result,
                    status="rejected",
                )
                try:
                    state = record_action(state, action)
                except ValueError as exc:
                    state = finish_session(
                        state,
                        status="blocked_by_policy",
                        stopping_reason=self._safe_error("invalid provider action", exc),
                    )
                    break
                state = finish_session(
                    state,
                    status="budget_exhausted",
                    stopping_reason=(
                        "provider response exceeded the session token or time budget"
                    ),
                )
                break
            if payload.action_name == "finish":
                action = self._action_record(
                    action_id,
                    payload,
                    normalized_arguments,
                    transport_result,
                    status="completed",
                )
                try:
                    state = record_action(state, action)
                    state = self._finish_from_payload(state, payload)
                except ValueError as exc:
                    state = finish_session(
                        state,
                        status="blocked_by_policy",
                        stopping_reason=self._safe_error("invalid final response", exc),
                    )
                break

            if payload.action_name == "propose_clarification":
                action = self._action_record(
                    action_id,
                    payload,
                    normalized_arguments,
                    transport_result,
                    status="completed",
                )
                try:
                    state = record_action(state, action)
                    question = self._question_from_payload(payload)
                    state = replace(
                        state,
                        usage=replace(
                            state.usage,
                            clarification_cycles=state.usage.clarification_cycles + 1,
                        ),
                    )
                    state = finish_session(
                        state,
                        status="needs_user_input",
                        stopping_reason=question.reason,
                        questions=(question,),
                    )
                except ValueError as exc:
                    state = finish_session(
                        state,
                        status="blocked_by_policy",
                        stopping_reason=self._safe_error("invalid clarification", exc),
                    )
                break

            try:
                capability_result = self._execute_action(
                    state,
                    payload.action_name,
                    normalized_arguments,
                )
                capability_evidence = self._collision_safe_evidence(
                    state,
                    capability_result.evidence,
                )
                new_ids = tuple(
                    item.evidence_id
                    for item in capability_evidence
                    if item.evidence_id not in state.allowed_evidence_ids
                )
                action = self._action_record(
                    action_id,
                    payload,
                    normalized_arguments,
                    transport_result,
                    status="completed",
                    result_ref=capability_result.result_ref,
                    added_evidence_ids=new_ids,
                )
                state = record_action(state, action)
                state = add_evidence(
                    state,
                    capability_evidence,
                    domain_result_ref=(
                        capability_result.result_ref
                        if capability_result.register_result_ref
                        else None
                    ),
                )
                last_result = capability_result.authoritative_result
                if payload.action_name == "retry_route_search":
                    state = replace(
                        state,
                        usage=replace(
                            state.usage,
                            search_expansions=state.usage.search_expansions + 1,
                        ),
                    )
                if payload.action_name == "apply_repair":
                    state = replace(
                        state,
                        usage=replace(
                            state.usage,
                            repair_attempts=state.usage.repair_attempts + 1,
                        ),
                    )
                no_progress_turns = 0 if new_ids else no_progress_turns + 1
            except ValueError as exc:
                state = finish_session(
                    state,
                    status="blocked_by_policy",
                    stopping_reason=self._safe_error("invalid capability action", exc),
                )
                break
            except Exception as exc:
                state = finish_session(
                    state,
                    status="domain_failed",
                    stopping_reason=self._safe_error("domain service failure", exc),
                )
                break

            if no_progress_turns >= self._policy.no_progress_turn_limit:
                state = finish_session(
                    state,
                    status="no_progress",
                    stopping_reason="two consecutive actions added no new evidence",
                )
                break

        return AssistanceRunResult(state=state, authoritative_result=last_result)

    def resume_with_condition_constraint(
        self,
        state: AssistanceSessionState,
        raw_value: Any,
    ) -> AssistanceRunResult:
        """Confirm one pending condition question and rerun canonically."""

        if state.status != "needs_user_input" or len(state.unresolved_questions) != 1:
            raise ValueError("session does not have exactly one pending clarification")
        question = state.unresolved_questions[0]
        if question.constraint_owner != "condition_registry":
            raise ValueError("pending clarification is not a condition constraint")
        resolution = normalize_condition_constraint(
            question.constraint_kind,  # type: ignore[arg-type]
            raw_value,
            provenance="confirmed_user",
        )
        if resolution.status != "resolved" or resolution.constraint is None:
            detail = ", ".join(resolution.warnings) or resolution.status
            raise ValueError(f"condition constraint was not uniquely resolved: {detail}")
        condition_constraints = state.request.condition_constraints.append(
            resolution.constraint
        )
        request = replace(
            state.request,
            condition_constraints=condition_constraints,
        )
        evidence_id = f"user.constraint.{resolution.constraint.constraint_id}"
        user_evidence = EvidenceItem(
            evidence_id=evidence_id,
            layer="user",
            source_id=question.question_id,
            payload_type="confirmed_condition_constraint",
            payload=resolution.constraint.to_dict(),
            provenance="user_confirmed",
            schema_versions={
                "condition_constraint": resolution.constraint.schema_version
            },
        )
        resumed = replace(
            state,
            request=request,
            unresolved_questions=(),
            status="active",
            stopping_reason=None,
        )
        resumed = add_evidence(resumed, (user_evidence,))
        return self.run(request, initial_state=resumed)

    def _allowed_actions(self, state: AssistanceSessionState) -> Tuple[str, ...]:
        declared = self._policy.actions_for(state.request.mode)
        if not state.domain_result_refs:
            initial = {
                "conditions": "recommend_conditions",
                "retro": "disconnect_target",
                "multistep": "plan_routes",
            }[state.request.mode]
            candidates = (initial, "propose_clarification", "finish")
            if state.request.mode in {"retro", "multistep"}:
                candidates = (
                    "audit_target",
                    *candidates,
                )
        else:
            candidates_by_mode = {
                "conditions": (
                    "inspect_condition_candidate",
                    "compare_condition_candidates",
                    "propose_clarification",
                    "finish",
                ),
                "retro": (
                    "audit_target",
                    "inspect_strategy",
                    "compare_strategies",
                    "inspect_strategy_conditions",
                    "apply_repair",
                    "verify_strategy",
                    "propose_clarification",
                    "finish",
                ),
                "multistep": (
                    "audit_target",
                    "inspect_route",
                    "inspect_route_step",
                    "search_step_precedents",
                    "compare_routes",
                    "apply_repair",
                    "verify_route",
                    "retry_route_search",
                    "propose_clarification",
                    "finish",
                ),
            }
            candidates = candidates_by_mode[state.request.mode]
        if (
            state.usage.clarification_cycles
            >= state.request.budget.max_clarification_cycles
        ):
            candidates = tuple(
                action for action in candidates if action != "propose_clarification"
            )
        if state.usage.search_expansions >= state.request.budget.max_search_expansions:
            candidates = tuple(
                action
                for action in candidates
                if action != "retry_route_search"
            )
        if state.usage.repair_attempts >= state.request.budget.max_repair_attempts:
            candidates = tuple(
                action for action in candidates if action != "apply_repair"
            )
        return tuple(action for action in candidates if action in declared)

    @staticmethod
    def _controller_packet(
        state: AssistanceSessionState,
        allowed_actions: Tuple[str, ...],
    ) -> Dict[str, Any]:
        return {
            "session_id": state.session_id,
            "mode": state.request.mode,
            "objective": state.request.objective,
            "structure_input": state.request.structure_input,
            "confirmed_constraints": [
                {
                    "constraint_id": item.constraint_id,
                    "owner": item.owner,
                    "kind": item.kind,
                    "value": item.value,
                    "provenance": item.provenance,
                    "hard": item.hard,
                }
                for item in state.request.constraints
                if item.provenance != "model_proposed"
            ],
            "allowed_actions": list(allowed_actions),
            "remaining_action_turns": state.remaining_action_turns,
            "domain_result_refs": list(state.domain_result_refs),
            "evidence": [
                {
                    "evidence_id": item.evidence_id,
                    "layer": item.layer,
                    "payload_type": item.payload_type,
                    "payload": dict(item.payload),
                    "provenance": item.provenance,
                    "uncertainty": item.uncertainty,
                }
                for item in state.evidence
            ],
        }

    def _normalized_arguments(
        self,
        payload: AssistanceActionPayload,
        state: AssistanceSessionState,
    ) -> Dict[str, Any]:
        arguments = payload.arguments.model_dump(exclude_none=True)
        if payload.action_name == "audit_target":
            if arguments:
                raise ValueError("audit_target does not accept model arguments")
        elif payload.action_name == "recommend_conditions":
            if arguments:
                raise ValueError("recommend_conditions does not accept model arguments")
            arguments["request_fingerprint"] = self._conditions.request_fingerprint(
                state.request
            )
        elif payload.action_name in {"disconnect_target", "plan_routes"}:
            if arguments:
                raise ValueError(f"{payload.action_name} does not accept model arguments")
            arguments["request_fingerprint"] = stable_assistance_id(
                "REQUEST",
                state.request,
            )
        elif payload.action_name in {
            "inspect_condition_candidate",
            "compare_condition_candidates",
            "inspect_strategy",
            "compare_strategies",
            "inspect_strategy_conditions",
            "inspect_route",
            "inspect_route_step",
            "search_step_precedents",
            "compare_routes",
            "apply_repair",
            "verify_strategy",
            "verify_route",
            "retry_route_search",
        }:
            requested_ref = arguments.get("result_ref")
            latest_ref = state.domain_result_refs[-1] if state.domain_result_refs else None
            if requested_ref is None:
                arguments["result_ref"] = latest_ref
            elif requested_ref != latest_ref:
                raise ValueError("action referenced a stale domain result")
        if payload.action_name == "apply_repair":
            allowed = {"result_ref", "proposal_id"}
            unknown = set(arguments) - allowed
            if unknown:
                raise ValueError(
                    f"apply_repair received unsupported arguments: {sorted(unknown)}"
                )
            proposal_id = arguments.get("proposal_id")
            if not isinstance(proposal_id, str) or not proposal_id:
                raise ValueError("apply_repair requires a proposal_id")
            cited = {
                item.evidence_id: item
                for item in state.evidence
                if item.evidence_id in payload.cited_evidence_ids
            }
            if not any(
                item.payload_type
                in {"route_repair_proposal", "single_step_repair_proposal"}
                and item.payload.get("status") == "actionable"
                and item.payload.get("proposal_id") == proposal_id
                for item in cited.values()
            ):
                raise ValueError(
                    "apply_repair must cite a matching actionable deterministic "
                    "repair proposal"
                )
        return arguments

    def _execute_action(
        self,
        state: AssistanceSessionState,
        action_name: str,
        arguments: Mapping[str, Any],
    ) -> CapabilityResult:
        if action_name == "audit_target":
            return self._chemistry.audit_target(state.request)
        if action_name == "recommend_conditions":
            return self._conditions.recommend_conditions(state.request)
        if action_name == "inspect_condition_candidate":
            alias = arguments.get("candidate_alias")
            if not isinstance(alias, str) or not alias:
                raise ValueError("candidate_alias is required")
            return self._conditions.inspect_condition_candidate(
                result_ref=str(arguments["result_ref"]),
                candidate_alias=alias,
            )
        if action_name == "compare_condition_candidates":
            aliases = arguments.get("candidate_aliases")
            if not isinstance(aliases, list):
                raise ValueError("candidate_aliases is required")
            return self._conditions.compare_condition_candidates(
                result_ref=str(arguments["result_ref"]),
                candidate_aliases=tuple(str(alias) for alias in aliases),
            )
        if action_name == "disconnect_target":
            return self._retro.disconnect_target(state.request)
        if action_name in {
            "inspect_strategy",
            "inspect_strategy_conditions",
        }:
            alias = arguments.get("strategy_alias")
            if not isinstance(alias, str) or not alias:
                raise ValueError("strategy_alias is required")
            if action_name == "inspect_strategy_conditions":
                return self._retro.inspect_strategy_conditions(
                    str(arguments["result_ref"]), alias
                )
            return self._retro.inspect_strategy(str(arguments["result_ref"]), alias)
        if action_name == "compare_strategies":
            aliases = arguments.get("strategy_aliases")
            if not isinstance(aliases, list):
                raise ValueError("strategy_aliases is required")
            return self._retro.compare_strategies(
                str(arguments["result_ref"]),
                tuple(str(alias) for alias in aliases),
            )
        if action_name == "plan_routes":
            return self._multistep.plan_routes(state.request)
        if action_name in {"inspect_route", "inspect_route_step"}:
            alias = arguments.get("route_alias")
            if not isinstance(alias, str) or not alias:
                raise ValueError("route_alias is required")
            step_index = arguments.get("step_index")
            if action_name == "inspect_route_step" and not isinstance(step_index, int):
                raise ValueError("step_index is required")
            return self._multistep.inspect_route(
                str(arguments["result_ref"]),
                alias,
                step_index=step_index if action_name == "inspect_route_step" else None,
            )
        if action_name == "search_step_precedents":
            alias = arguments.get("route_alias")
            step_index = arguments.get("step_index")
            if not isinstance(alias, str) or not alias:
                raise ValueError("route_alias is required")
            if not isinstance(step_index, int):
                raise ValueError("step_index is required")
            return self._multistep.search_step_precedents(
                str(arguments["result_ref"]),
                alias,
                step_index=step_index,
            )
        if action_name == "compare_routes":
            aliases = arguments.get("route_aliases")
            if not isinstance(aliases, list):
                raise ValueError("route_aliases is required")
            return self._multistep.compare_routes(
                str(arguments["result_ref"]),
                tuple(str(alias) for alias in aliases),
            )
        if action_name == "apply_repair":
            proposal_id = arguments.get("proposal_id")
            if not isinstance(proposal_id, str) or not proposal_id:
                raise ValueError("proposal_id is required")
            if state.request.mode == "retro":
                return self._retro.apply_repair(
                    state.request,
                    str(arguments["result_ref"]),
                    proposal_id=proposal_id,
                )
            return self._multistep.apply_repair(
                state.request,
                str(arguments["result_ref"]),
                proposal_id=proposal_id,
            )
        if action_name == "verify_strategy":
            alias = arguments.get("strategy_alias")
            if not isinstance(alias, str) or not alias:
                raise ValueError("strategy_alias is required")
            return self._retro.verify_strategy(
                str(arguments["result_ref"]), alias
            )
        if action_name == "verify_route":
            alias = arguments.get("route_alias")
            if not isinstance(alias, str) or not alias:
                raise ValueError("route_alias is required")
            return self._multistep.verify_route(
                str(arguments["result_ref"]),
                alias,
            )
        if action_name == "retry_route_search":
            delta = RouteSearchPolicyDelta(
                max_depth_delta=int(arguments.get("search_depth_delta", 0)),
                beam_width_delta=int(arguments.get("beam_width_delta", 0)),
                max_expansions_delta=int(arguments.get("max_expansions_delta", 0)),
            )
            return self._multistep.retry_route_search(
                state.request,
                str(arguments["result_ref"]),
                delta,
            )
        raise ValueError(f"unsupported assistance capability: {action_name}")

    @staticmethod
    def _action_record(
        action_id: str,
        payload: AssistanceActionPayload,
        normalized_arguments: Mapping[str, Any],
        result: AssistanceTransportResult,
        *,
        status: str,
        result_ref: str | None = None,
        added_evidence_ids: Tuple[str, ...] = (),
    ) -> ActionRecord:
        return ActionRecord(
            action_id=action_id,
            action_name=payload.action_name,
            normalized_arguments=normalized_arguments,
            cited_evidence_ids=tuple(payload.cited_evidence_ids),
            status=status,
            result_ref=result_ref,
            added_evidence_ids=added_evidence_ids,
            decision_summary=payload.decision_summary,
            provider_attempts=result.attempts,
            input_tokens=result.input_tokens,
            output_tokens=result.output_tokens,
            elapsed_ms=result.elapsed_ms,
        )

    @staticmethod
    def _question_from_payload(
        payload: AssistanceActionPayload,
    ) -> ClarificationQuestion:
        values = payload.arguments
        required = {
            "question_id": values.question_id,
            "prompt": values.prompt,
            "constraint_owner": values.constraint_owner,
            "constraint_kind": values.constraint_kind,
            "reason": values.reason,
        }
        missing = [name for name, value in required.items() if not value]
        if missing:
            raise ValueError(f"clarification fields are required: {sorted(missing)!r}")
        return ClarificationQuestion(
            question_id=str(values.question_id),
            prompt=str(values.prompt),
            constraint_owner=str(values.constraint_owner),
            constraint_kind=str(values.constraint_kind),
            reason=str(values.reason),
            evidence_ids=tuple(payload.cited_evidence_ids),
        )

    @staticmethod
    def _finish_from_payload(
        state: AssistanceSessionState,
        payload: AssistanceActionPayload,
    ) -> AssistanceSessionState:
        terminal_status = payload.arguments.terminal_status
        reason = payload.arguments.stopping_reason
        if terminal_status is None or not reason:
            raise ValueError("finish requires terminal_status and stopping_reason")
        claims = tuple(
            AdvisoryClaim(
                claim_type=item.claim_type,
                subject_id=item.subject_id,
                severity=item.severity,
                statement=item.statement,
                evidence_ids=tuple(item.evidence_ids),
                uncertainty=item.uncertainty,
                suggested_user_action=item.suggested_user_action,
            )
            for item in payload.claims
        )
        if terminal_status == "completed" and not claims:
            raise ValueError("completed assistance requires at least one claim")
        return finish_session(
            state,
            status=terminal_status,
            stopping_reason=reason,
            claims=claims,
        )

    def _last_result(self, state: AssistanceSessionState) -> object | None:
        capabilities = self._capabilities_for(state.request.mode)
        if not state.domain_result_refs or capabilities is None:
            return None
        try:
            return capabilities.result(state.domain_result_refs[-1])
        except ValueError:
            return None

    def _capabilities_for(self, mode: str) -> object | None:
        return {
            "conditions": self._conditions,
            "retro": self._retro,
            "multistep": self._multistep,
        }.get(mode)

    @staticmethod
    def _collision_safe_evidence(
        state: AssistanceSessionState,
        evidence: Tuple[EvidenceItem, ...],
    ) -> Tuple[EvidenceItem, ...]:
        """Namespace changed evidence while retaining the immutable prior item."""

        existing = {item.evidence_id: item for item in state.evidence}
        normalized = []
        for item in evidence:
            prior = existing.get(item.evidence_id)
            if prior is None or prior == item:
                normalized.append(item)
                continue
            source_token = item.source_id.replace(":", "-")[-20:]
            candidate_id = f"revision-{source_token}.{item.evidence_id}"
            suffix = 2
            while candidate_id in existing and existing[candidate_id] != item:
                candidate_id = (
                    f"revision-{source_token}-{suffix}.{item.evidence_id}"
                )
                suffix += 1
            revised = replace(item, evidence_id=candidate_id)
            normalized.append(revised)
            existing[candidate_id] = revised
        return tuple(normalized)

    @staticmethod
    def _budget_reason(state: AssistanceSessionState) -> str | None:
        budget = state.request.budget
        if state.usage.action_turns >= budget.max_action_turns:
            return "maximum action turns exhausted"
        if state.usage.total_tokens >= budget.max_total_tokens:
            return "maximum provider tokens exhausted"
        if state.usage.elapsed_ms >= budget.max_elapsed_ms:
            return "maximum elapsed provider time exhausted"
        return None

    @staticmethod
    def _safe_error(prefix: str, exc: Exception) -> str:
        message = str(exc).replace(chr(10), " ")[:300]
        return f"{prefix}: {_SECRET_PATTERN.sub('[REDACTED]', message)}"
