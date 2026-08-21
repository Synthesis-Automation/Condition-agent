"""Deterministic rendering for advisory assistance state."""

from __future__ import annotations

from .contracts import AssistanceSessionState


def render_assistance(state: AssistanceSessionState) -> str:
    """Render claims, uncertainty, questions, and the stopping reason."""

    lines = [f"Assistance status: {state.status}."]
    if state.stopping_reason:
        lines.append(state.stopping_reason)
    for claim in state.claims:
        lines.append(
            f"- [{claim.severity}] {claim.statement} "
            f"(evidence: {', '.join(claim.evidence_ids)})"
        )
        if claim.uncertainty:
            lines.append(f"  Uncertainty: {claim.uncertainty}")
    for question in state.unresolved_questions:
        lines.append(f"Question: {question.prompt}")
    return "\n".join(lines)
