"""Agent-level validation gate for final reaction decisions."""

from __future__ import annotations

from typing import Any, Dict, List

from .contracts import ValidationReport


class AgentDecisionValidator:
    """Validate analysis payloads before finalizing agent output."""

    def __init__(self, *, min_confidence: float = 0.5) -> None:
        self.min_confidence = min_confidence

    def validate_analysis(self, analysis: Dict[str, Any]) -> ValidationReport:
        """Validate the analysis and return a final decision."""
        issues: List[str] = []
        decision = dict(analysis.get("decision") or {})
        candidate_rows = analysis.get("taxonomy_candidates") or []
        candidate_ids = {
            str(row.get("reaction_type"))
            for row in candidate_rows
            if isinstance(row, dict) and row.get("reaction_type")
        }

        reaction_type = str(decision.get("reaction_type") or "unknown")
        confidence = float(decision.get("confidence") or 0.0)

        if reaction_type == "unknown":
            return ValidationReport(
                passed=True,
                issues=[],
                final_decision={
                    "reaction_type": "unknown",
                    "confidence": confidence,
                    "source": decision.get("source", "validator"),
                    "rationale": decision.get("rationale", "No confident taxonomy candidate."),
                },
            )

        if reaction_type not in candidate_ids:
            issues.append("decision_not_in_candidate_list")
        if confidence < self.min_confidence:
            issues.append("decision_confidence_below_threshold")

        if issues:
            return ValidationReport(
                passed=False,
                issues=issues,
                final_decision={
                    "reaction_type": "unknown",
                    "confidence": 0.0,
                    "source": "validator_fallback",
                    "rationale": "Validation failed: " + "; ".join(issues),
                },
            )

        return ValidationReport(
            passed=True,
            issues=[],
            final_decision={
                "reaction_type": reaction_type,
                "confidence": confidence,
                "source": decision.get("source", "deterministic"),
                "rationale": decision.get("rationale", "Validated deterministic decision."),
            },
        )
