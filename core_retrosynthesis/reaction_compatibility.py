"""Planner policy for graph-derived reaction-regime compatibility."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Tuple

from reactive_taxonomy import (
    ReactionCompatibilityAssessment,
    assess_reaction_compatibility,
)


REACTION_COMPATIBILITY_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_compatibility_policy.v1.json"
)
_SEVERITIES = ("low", "medium", "high", "critical")
_DISPOSITIONS = ("pass", "warn", "demote", "reject")


@dataclass(frozen=True)
class ReactionCompatibilityPolicy:
    """Validated mapping from reaction-regime severity to planner action."""

    definition_id: str
    schema_version: str
    severity_actions: Tuple[Tuple[str, str, int], ...]
    disposition_order: Tuple[str, ...]

    def action(self, severity: str) -> tuple[str, int]:
        """Return disposition and structural-band penalty for a severity."""

        actions = {
            key: (disposition, penalty)
            for key, disposition, penalty in self.severity_actions
        }
        try:
            return actions[severity]
        except KeyError as error:
            raise ValueError(
                f"unsupported reaction compatibility severity: {severity}"
            ) from error


@dataclass(frozen=True)
class ReactionCompatibilityResult:
    """Aggregated reaction-regime compatibility result for one candidate."""

    assessments: Tuple[ReactionCompatibilityAssessment, ...]
    disposition: str
    structural_band_penalty: int
    warning_strength: str | None
    policy_definition_id: str
    policy_schema_version: str

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible result."""

        return asdict(self)


def validate_reaction_compatibility_policy(payload: Mapping[str, Any]) -> None:
    """Validate the planner reaction-compatibility policy."""

    if payload.get("definition_id") != "reaction_compatibility_policy.v1":
        raise ValueError("unexpected reaction compatibility policy ID")
    if payload.get("schema_version") != "1.0":
        raise ValueError("unsupported reaction compatibility policy schema")
    actions = payload.get("severity_actions")
    if not isinstance(actions, Mapping) or set(actions) != set(_SEVERITIES):
        raise ValueError("reaction compatibility severity actions are incomplete")
    for severity in _SEVERITIES:
        action = actions[severity]
        if not isinstance(action, Mapping):
            raise ValueError(
                f"reaction compatibility action {severity} must be an object"
            )
        if action.get("disposition") not in _DISPOSITIONS:
            raise ValueError(
                f"invalid reaction compatibility disposition for {severity}"
            )
        penalty = action.get("structural_band_penalty")
        if isinstance(penalty, bool) or not isinstance(penalty, int) or penalty < 0:
            raise ValueError(
                f"invalid reaction compatibility band penalty for {severity}"
            )
    aggregation = payload.get("aggregation")
    if not isinstance(aggregation, Mapping):
        raise ValueError("reaction compatibility aggregation must be an object")
    if aggregation.get("band_penalty") != "maximum":
        raise ValueError("reaction compatibility penalties must aggregate by maximum")
    if tuple(aggregation.get("disposition_order") or ()) != _DISPOSITIONS:
        raise ValueError("reaction compatibility disposition order is unsupported")


@lru_cache(maxsize=1)
def load_reaction_compatibility_policy() -> ReactionCompatibilityPolicy:
    """Load the canonical reaction-regime planner policy."""

    payload = json.loads(
        REACTION_COMPATIBILITY_POLICY_PATH.read_text(encoding="utf-8")
    )
    validate_reaction_compatibility_policy(payload)
    actions = payload["severity_actions"]
    return ReactionCompatibilityPolicy(
        definition_id=str(payload["definition_id"]),
        schema_version=str(payload["schema_version"]),
        severity_actions=tuple(
            (
                severity,
                str(actions[severity]["disposition"]),
                int(actions[severity]["structural_band_penalty"]),
            )
            for severity in _SEVERITIES
        ),
        disposition_order=_DISPOSITIONS,
    )


@lru_cache(maxsize=20_000)
def assess_candidate_reaction_compatibility(
    reaction_smiles: str,
) -> ReactionCompatibilityResult:
    """Assess one candidate and aggregate definition-derived planner actions."""

    policy = load_reaction_compatibility_policy()
    assessments = assess_reaction_compatibility(reaction_smiles)
    if not assessments:
        return ReactionCompatibilityResult(
            assessments=(),
            disposition="pass",
            structural_band_penalty=0,
            warning_strength=None,
            policy_definition_id=policy.definition_id,
            policy_schema_version=policy.schema_version,
        )
    actions = tuple(
        policy.action(assessment.intrinsic_severity)
        for assessment in assessments
    )
    disposition = max(
        (action[0] for action in actions),
        key=policy.disposition_order.index,
    )
    warning_strength = (
        "strong"
        if any(assessment.warning_strength == "strong" for assessment in assessments)
        else "advisory"
    )
    return ReactionCompatibilityResult(
        assessments=assessments,
        disposition=disposition,
        structural_band_penalty=max(action[1] for action in actions),
        warning_strength=warning_strength,
        policy_definition_id=policy.definition_id,
        policy_schema_version=policy.schema_version,
    )


__all__ = [
    "REACTION_COMPATIBILITY_POLICY_PATH",
    "ReactionCompatibilityPolicy",
    "ReactionCompatibilityResult",
    "assess_candidate_reaction_compatibility",
    "load_reaction_compatibility_policy",
    "validate_reaction_compatibility_policy",
]
