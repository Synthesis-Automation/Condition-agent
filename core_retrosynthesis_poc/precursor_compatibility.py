"""Planner policy for transparent intramolecular compatibility warnings."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Tuple

from reactive_taxonomy.reactive_pair_interactions import (
    ReactivePairInteractionAssessment,
    assess_reactive_pair_interactions,
)


PRECURSOR_COMPATIBILITY_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "precursor_compatibility_policy.v1.json"
)
_SEVERITIES = ("low", "medium", "high", "critical")
_DISPOSITIONS = ("pass", "warn", "demote", "reject")


@dataclass(frozen=True)
class PrecursorCompatibilityPolicy:
    """Validated immutable mapping from chemistry severity to planner action."""

    definition_id: str
    schema_version: str
    severity_actions: Tuple[Tuple[str, str, int], ...]
    disposition_order: Tuple[str, ...]

    def action(self, severity: str) -> tuple[str, int]:
        """Return disposition and band penalty for a chemistry severity."""

        actions = {
            key: (disposition, penalty)
            for key, disposition, penalty in self.severity_actions
        }
        try:
            return actions[severity]
        except KeyError as error:
            raise ValueError(f"unsupported compatibility severity: {severity}") from error


@dataclass(frozen=True)
class PrecursorCompatibilityResult:
    """Aggregated, auditable compatibility result for precursor components."""

    assessments: Tuple[ReactivePairInteractionAssessment, ...]
    disposition: str
    structural_band_penalty: int
    warning_strength: str | None
    policy_definition_id: str
    policy_schema_version: str

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible result."""

        return asdict(self)


def validate_precursor_compatibility_policy(payload: Mapping[str, Any]) -> None:
    """Validate a planner compatibility policy without executing chemistry."""

    if payload.get("definition_id") != "precursor_compatibility_policy.v1":
        raise ValueError("unexpected precursor compatibility policy ID")
    if payload.get("schema_version") != "1.0":
        raise ValueError("unsupported precursor compatibility policy schema")
    scope = payload.get("scope_policy")
    if not isinstance(scope, Mapping) or dict(scope) != {
        "same_component": "apply",
        "different_components": "no_intrinsic_molecular_penalty",
    }:
        raise ValueError("precursor compatibility policy must be component scoped")
    actions = payload.get("severity_actions")
    if not isinstance(actions, Mapping) or set(actions) != set(_SEVERITIES):
        raise ValueError("precursor compatibility severity actions are incomplete")
    for severity in _SEVERITIES:
        action = actions[severity]
        if not isinstance(action, Mapping):
            raise ValueError(f"compatibility action {severity} must be an object")
        if action.get("disposition") not in _DISPOSITIONS:
            raise ValueError(f"invalid compatibility disposition for {severity}")
        penalty = action.get("structural_band_penalty")
        if isinstance(penalty, bool) or not isinstance(penalty, int) or penalty < 0:
            raise ValueError(f"invalid compatibility band penalty for {severity}")
    aggregation = payload.get("aggregation")
    if not isinstance(aggregation, Mapping):
        raise ValueError("precursor compatibility aggregation must be an object")
    if aggregation.get("band_penalty") != "maximum":
        raise ValueError("compatibility penalties must aggregate by maximum")
    if tuple(aggregation.get("disposition_order") or ()) != _DISPOSITIONS:
        raise ValueError("compatibility disposition order is unsupported")


@lru_cache(maxsize=1)
def load_precursor_compatibility_policy() -> PrecursorCompatibilityPolicy:
    """Load the canonical planner action policy."""

    payload = json.loads(
        PRECURSOR_COMPATIBILITY_POLICY_PATH.read_text(encoding="utf-8")
    )
    validate_precursor_compatibility_policy(payload)
    actions = payload["severity_actions"]
    return PrecursorCompatibilityPolicy(
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
def assess_precursor_compatibility(
    precursor_smiles: str,
) -> PrecursorCompatibilityResult:
    """Assess only incompatible pairs contained in the same precursor component."""

    policy = load_precursor_compatibility_policy()
    assessments = assess_reactive_pair_interactions(precursor_smiles)
    if not assessments:
        return PrecursorCompatibilityResult(
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
    return PrecursorCompatibilityResult(
        assessments=assessments,
        disposition=disposition,
        structural_band_penalty=max(action[1] for action in actions),
        warning_strength=warning_strength,
        policy_definition_id=policy.definition_id,
        policy_schema_version=policy.schema_version,
    )


__all__ = [
    "PrecursorCompatibilityPolicy",
    "PrecursorCompatibilityResult",
    "assess_precursor_compatibility",
    "load_precursor_compatibility_policy",
    "validate_precursor_compatibility_policy",
]
