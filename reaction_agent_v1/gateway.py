"""Gateway orchestration for reaction agent sessions."""

from __future__ import annotations

from typing import Any, Dict, Optional

from .contracts import AgentRunResult, SessionState
from .evidence import ReactionEvidence
from .memory import SessionMemory
from .planner import DynamicPlanner
from .tool_registry import ToolRegistry, build_default_registry


class ReactionAgentGateway:
    """Session-aware gateway for planner-tool execution loops."""

    def __init__(
        self,
        *,
        memory: Optional[SessionMemory] = None,
        planner: Optional[DynamicPlanner] = None,
        registry: Optional[ToolRegistry] = None,
        min_confidence: float = 0.5,
    ) -> None:
        self.memory = memory or SessionMemory()
        self.planner = planner or DynamicPlanner()
        self.registry = registry or build_default_registry(min_confidence=min_confidence)

    def run(
        self,
        reaction_smiles: str,
        *,
        session_id: Optional[str] = None,
        max_steps: int = 12,
    ) -> AgentRunResult:
        """Run one analysis loop for a reaction."""
        state = self.memory.create_session(reaction_smiles=reaction_smiles, session_id=session_id)
        state.context["evidence"] = ReactionEvidence(reaction_smiles=reaction_smiles).to_dict()
        self._execute_loop(state=state, max_steps=max_steps)
        return self._build_result(state)

    def _execute_loop(self, *, state: SessionState, max_steps: int) -> None:
        for step_index in range(max_steps):
            plan = self.planner.next_step(state)
            if plan.action == "stop":
                self.memory.append_trace(
                    state,
                    step_index=step_index,
                    action=plan.action,
                    tool_name=plan.tool_name,
                    status="ok",
                    message=plan.reason,
                )
                break

            if plan.action == "finalize":
                self._finalize_state(state)
                self.memory.append_trace(
                    state,
                    step_index=step_index,
                    action=plan.action,
                    tool_name=plan.tool_name,
                    status="ok",
                    message=plan.reason,
                )
                continue

            tool_result = self.registry.execute(plan.tool_name or "", **self._tool_kwargs(plan.action, state))
            if not tool_result.ok:
                state.status = "failed"
                self.memory.append_trace(
                    state,
                    step_index=step_index,
                    action=plan.action,
                    tool_name=plan.tool_name,
                    status="error",
                    message=tool_result.error or "tool failed",
                )
                break

            self._persist_tool_payload(state=state, action=plan.action, payload=tool_result.payload)
            self.memory.append_trace(
                state,
                step_index=step_index,
                action=plan.action,
                tool_name=plan.tool_name,
                status="ok",
                message=f"{plan.reason} (latency_ms={tool_result.latency_ms:.1f})",
            )

        if state.status == "running":
            self._finalize_state(state)

    def _tool_kwargs(self, action: str, state: SessionState) -> Dict[str, Any]:
        evidence = self._load_evidence(state)
        if action == "reaction_diff":
            return {"reaction_smiles": state.reaction_smiles}
        if action == "fallback_candidates":
            return {
                "reaction_smiles": state.reaction_smiles,
                "evidence": evidence.to_dict(),
            }
        if action == "validate":
            return {"evidence": evidence.to_dict()}
        if action == "confidence_calibration":
            return {
                "evidence": evidence.to_dict(),
                "validation": dict(evidence.validation),
            }
        if action == "llm_rerank":
            return {
                "evidence": evidence.to_dict(),
                "validation": dict(evidence.validation),
            }
        if action == "precedent_lookup":
            return {
                "reaction_smiles": state.reaction_smiles,
                "evidence": evidence.to_dict(),
            }
        if action == "coverage":
            return {
                "reaction_smiles": state.reaction_smiles,
                "evidence": evidence.to_dict(),
            }
        return {}

    def _persist_tool_payload(self, *, state: SessionState, action: str, payload: Dict[str, Any]) -> None:
        evidence = self._load_evidence(state)
        if action == "reaction_diff":
            evidence.merge_diff_payload(payload)
            self._save_evidence(state, evidence)
            return
        if action == "fallback_candidates":
            evidence.tool_artifacts["fallback_candidates"] = payload
            recovered = list(payload.get("recovered_candidates") or [])
            if recovered:
                evidence.merge_fallback_candidates(recovered)
            self._save_evidence(state, evidence)
            return
        if action == "validate":
            evidence.validation = dict(payload)
            validation_final = payload.get("final_decision") or {}
            if validation_final:
                evidence.final_decision = dict(validation_final)
            self._save_evidence(state, evidence)
            return
        if action == "confidence_calibration":
            evidence.tool_artifacts["confidence_calibration"] = payload
            self._save_evidence(state, evidence)
            return
        if action == "llm_rerank":
            evidence.tool_artifacts["llm_rerank"] = payload
            self._save_evidence(state, evidence)
            return
        if action == "precedent_lookup":
            evidence.tool_artifacts["precedent_lookup"] = payload
            self._save_evidence(state, evidence)
            return
        if action == "coverage":
            evidence.coverage_suggestions = list(payload.get("suggestions") or [])
            self._save_evidence(state, evidence)

    def _finalize_state(self, state: SessionState) -> None:
        evidence = self._load_evidence(state)
        if not evidence.final_decision:
            validation_final = (evidence.validation or {}).get("final_decision") or {}
            if validation_final:
                evidence.final_decision = dict(validation_final)
            else:
                evidence.final_decision = dict(
                    evidence.provisional_decision
                    or {"reaction_type": "unknown", "confidence": 0.0, "source": "finalize_fallback"}
                )
        self._save_evidence(state, evidence)
        state.status = "completed"

    def _build_result(self, state: SessionState) -> AgentRunResult:
        evidence = self._load_evidence(state)
        return AgentRunResult(
            session_id=state.session_id,
            reaction_smiles=state.reaction_smiles,
            status=state.status,
            final_decision=dict(evidence.final_decision or evidence.provisional_decision),
            evidence=evidence.to_dict(),
            analysis=evidence.to_analysis_view(),
            validation=dict(evidence.validation),
            tool_artifacts=dict(evidence.tool_artifacts),
            coverage_suggestions=list(evidence.coverage_suggestions),
            trace=[event for event in state.trace],
        )

    def export_session(self, session_id: str) -> Dict[str, Any]:
        """Export session state as dictionary."""
        return self.memory.export_session(session_id)

    def _load_evidence(self, state: SessionState) -> ReactionEvidence:
        payload = state.context.get("evidence")
        if isinstance(payload, dict) and payload.get("reaction_smiles"):
            return ReactionEvidence.from_dict(payload)
        evidence = ReactionEvidence(reaction_smiles=state.reaction_smiles)
        self._save_evidence(state, evidence)
        return evidence

    def _save_evidence(self, state: SessionState, evidence: ReactionEvidence) -> None:
        state.context["evidence"] = evidence.to_dict()
