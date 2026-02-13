"""Dynamic planner for reaction agent tool orchestration."""

from __future__ import annotations

from .contracts import PlanStep, SessionState


class DynamicPlanner:
    """Simple deterministic planner for foundation workflows."""

    def next_step(self, state: SessionState) -> PlanStep:
        """Decide the next action from current session context."""
        if "analysis" not in state.context:
            return PlanStep(
                action="analyze",
                tool_name="analyze_reaction",
                reason="Need deterministic taxonomy analysis before any decision.",
            )

        if "validation" not in state.context:
            return PlanStep(
                action="validate",
                tool_name="validate_decision",
                reason="Validate the selected decision against candidate and confidence gates.",
            )

        final_decision = state.context.get("final_decision") or {}
        if not final_decision:
            validation = state.context.get("validation") or {}
            validated_final = validation.get("final_decision") or {}
            if validated_final:
                state.context["final_decision"] = validated_final
                final_decision = validated_final

        if final_decision.get("reaction_type") == "unknown" and "coverage_suggestions" not in state.context:
            return PlanStep(
                action="coverage",
                tool_name="coverage_advice",
                reason="Unknown decision should trigger taxonomy/tool coverage expansion advice.",
            )

        if state.status != "completed":
            return PlanStep(
                action="finalize",
                tool_name=None,
                reason="All required stages are complete; finalize the run.",
            )

        return PlanStep(action="stop", tool_name=None, reason="Session already completed.")
