"""Gateway orchestration for reaction agent sessions."""

from __future__ import annotations

from typing import Any, Dict, Optional

from .contracts import AgentRunResult, SessionState
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
        if action == "analyze":
            return {"reaction_smiles": state.reaction_smiles}
        if action == "fallback_candidates":
            return {
                "reaction_smiles": state.reaction_smiles,
                "analysis": state.context.get("analysis", {}),
            }
        if action == "validate":
            return {"analysis": state.context.get("analysis", {})}
        if action == "confidence_calibration":
            return {
                "analysis": state.context.get("analysis", {}),
                "validation": state.context.get("validation", {}),
            }
        if action == "llm_rerank":
            return {
                "analysis": state.context.get("analysis", {}),
                "validation": state.context.get("validation", {}),
            }
        if action == "precedent_lookup":
            return {
                "reaction_smiles": state.reaction_smiles,
                "analysis": state.context.get("analysis", {}),
            }
        if action == "coverage":
            return {
                "reaction_smiles": state.reaction_smiles,
                "analysis": state.context.get("analysis", {}),
            }
        return {}

    def _persist_tool_payload(self, *, state: SessionState, action: str, payload: Dict[str, Any]) -> None:
        artifacts = state.context.setdefault("tool_artifacts", {})
        if action == "analyze":
            state.context["analysis"] = payload
            return
        if action == "fallback_candidates":
            state.context["fallback_candidates"] = payload
            artifacts["fallback_candidates"] = payload
            return
        if action == "validate":
            state.context["validation"] = payload
            validation_final = payload.get("final_decision") or {}
            if validation_final:
                state.context["final_decision"] = validation_final
            return
        if action == "confidence_calibration":
            state.context["confidence_calibration"] = payload
            artifacts["confidence_calibration"] = payload
            return
        if action == "llm_rerank":
            state.context["llm_rerank"] = payload
            artifacts["llm_rerank"] = payload
            return
        if action == "precedent_lookup":
            state.context["precedent_lookup"] = payload
            artifacts["precedent_lookup"] = payload
            return
        if action == "coverage":
            state.context["coverage_suggestions"] = payload.get("suggestions", [])

    def _finalize_state(self, state: SessionState) -> None:
        if "final_decision" not in state.context:
            validation = state.context.get("validation") or {}
            validation_final = validation.get("final_decision") or {}
            if validation_final:
                state.context["final_decision"] = validation_final
            else:
                analysis = state.context.get("analysis") or {}
                state.context["final_decision"] = analysis.get("decision", {"reaction_type": "unknown", "confidence": 0.0})
        state.status = "completed"

    def _build_result(self, state: SessionState) -> AgentRunResult:
        return AgentRunResult(
            session_id=state.session_id,
            reaction_smiles=state.reaction_smiles,
            status=state.status,
            final_decision=dict(state.context.get("final_decision") or {}),
            analysis=dict(state.context.get("analysis") or {}),
            validation=dict(state.context.get("validation") or {}),
            tool_artifacts=dict(state.context.get("tool_artifacts") or {}),
            coverage_suggestions=list(state.context.get("coverage_suggestions") or []),
            trace=[event for event in state.trace],
        )

    def export_session(self, session_id: str) -> Dict[str, Any]:
        """Export session state as dictionary."""
        return self.memory.export_session(session_id)
