"""Dynamic planner for reaction agent tool orchestration."""

from __future__ import annotations

from .contracts import PlanStep, SessionState
from .evidence import ReactionEvidence


def _load_evidence(state: SessionState) -> ReactionEvidence:
    payload = state.context.get("evidence")
    if isinstance(payload, dict) and payload.get("reaction_smiles"):
        return ReactionEvidence.from_dict(payload)
    evidence = ReactionEvidence(reaction_smiles=state.reaction_smiles)
    state.context["evidence"] = evidence.to_dict()
    return evidence


def _save_evidence(state: SessionState, evidence: ReactionEvidence) -> None:
    state.context["evidence"] = evidence.to_dict()


class DynamicPlanner:
    """Deterministic planner driven by evidence completeness."""

    def next_step(self, state: SessionState) -> PlanStep:
        """Decide the next action from current session context."""
        evidence = _load_evidence(state)

        if not evidence.has_diff():
            return PlanStep(
                action="reaction_diff",
                tool_name="reaction_diff",
                reason="Need canonical reaction diff evidence before classification and validation.",
            )

        if not evidence.has_validation():
            if not evidence.taxonomy_candidates and "fallback_candidates" not in evidence.tool_artifacts:
                return PlanStep(
                    action="fallback_candidates",
                    tool_name="fallback_candidate_retrieval",
                    reason="No deterministic candidates; query fallback candidate retrieval.",
                )
            if evidence.taxonomy_candidates and "fallback_candidates" not in evidence.tool_artifacts:
                evidence.tool_artifacts["fallback_candidates"] = {
                    "status": "skipped",
                    "reason": "deterministic_candidates_available",
                }
                _save_evidence(state, evidence)
            return PlanStep(
                action="validate",
                tool_name="validate_decision",
                reason="Validate provisional decision against evidence and candidate constraints.",
            )

        if "confidence_calibration" not in evidence.tool_artifacts:
            return PlanStep(
                action="confidence_calibration",
                tool_name="confidence_calibrator",
                reason="Run confidence calibration stage (placeholder contract in foundation mode).",
            )

        if "llm_rerank" not in evidence.tool_artifacts:
            return PlanStep(
                action="llm_rerank",
                tool_name="llm_rerank_constrained",
                reason="Run constrained LLM rerank stage (placeholder contract in foundation mode).",
            )

        if "precedent_lookup" not in evidence.tool_artifacts:
            return PlanStep(
                action="precedent_lookup",
                tool_name="precedent_lookup",
                reason="Run precedent lookup stage (placeholder contract in foundation mode).",
            )

        if evidence.decision_reaction_type() == "unknown" and not evidence.coverage_suggestions:
            return PlanStep(
                action="coverage",
                tool_name="coverage_advice",
                reason="Unknown decision should trigger taxonomy/tool coverage expansion advice.",
            )

        if state.status != "completed":
            return PlanStep(
                action="finalize",
                tool_name=None,
                reason="All required evidence stages are complete; finalize the run.",
            )

        return PlanStep(action="stop", tool_name=None, reason="Session already completed.")
