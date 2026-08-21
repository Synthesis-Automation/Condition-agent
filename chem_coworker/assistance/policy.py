"""Validated loading for versioned assistance policy definitions."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Mapping, Tuple

from .contracts import AssistanceBudget, AssistanceMode


_DEFINITIONS_DIR = Path(__file__).resolve().parent.parent / "definitions"
_DEFAULT_POLICY = _DEFINITIONS_DIR / "assistance_policy.v1.json"


@dataclass(frozen=True)
class AssistancePolicy:
    """Closed capability surface and validated controller bounds."""

    definition_id: str
    definition_version: str
    schema_version: str
    default_budget: AssistanceBudget
    maximum_budget: AssistanceBudget
    allowed_actions: Mapping[str, Tuple[str, ...]]
    no_progress_turn_limit: int
    default_rollout_state: str
    rollout_states: Tuple[str, ...]

    def actions_for(self, mode: AssistanceMode) -> Tuple[str, ...]:
        """Return statically declared actions for one assistance mode."""

        return self.allowed_actions[mode]

    def validate_budget(self, budget: AssistanceBudget) -> None:
        """Reject any caller budget that exceeds the versioned maxima."""

        requested = budget.__dict__
        maximum = self.maximum_budget.__dict__
        exceeded = [name for name, value in requested.items() if value > maximum[name]]
        if exceeded:
            raise ValueError(f"assistance budget exceeds policy: {sorted(exceeded)!r}")


def _budget(value: object, *, field_name: str) -> AssistanceBudget:
    if not isinstance(value, dict):
        raise ValueError(f"{field_name} must be an object")
    try:
        return AssistanceBudget(**value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field_name}: {exc}") from exc


@lru_cache(maxsize=4)
def load_assistance_policy(path: str | Path = _DEFAULT_POLICY) -> AssistancePolicy:
    """Load and strictly validate one declarative assistance policy."""

    definition_path = Path(path)
    data = json.loads(definition_path.read_text(encoding="utf-8"))
    required = {
        "definition_id",
        "definition_version",
        "schema_version",
        "default_budget",
        "maximum_budget",
        "allowed_actions",
        "no_progress_turn_limit",
        "default_rollout_state",
        "rollout_states",
    }
    if set(data) != required:
        raise ValueError(
            f"assistance policy fields must be exactly {sorted(required)!r}"
        )
    actions = data["allowed_actions"]
    if not isinstance(actions, dict) or set(actions) != {
        "conditions",
        "retro",
        "multistep",
    }:
        raise ValueError("allowed_actions must define all assistance modes")
    normalized: Dict[str, Tuple[str, ...]] = {}
    for mode, values in actions.items():
        if not isinstance(values, list) or not values:
            raise ValueError(f"allowed actions for {mode} must be a non-empty list")
        if any(not isinstance(value, str) or not value for value in values):
            raise ValueError(f"allowed actions for {mode} must be named strings")
        if len(values) != len(set(values)):
            raise ValueError(f"allowed actions for {mode} must be unique")
        normalized[mode] = tuple(values)
    default_budget = _budget(data["default_budget"], field_name="default_budget")
    maximum_budget = _budget(data["maximum_budget"], field_name="maximum_budget")
    policy = AssistancePolicy(
        definition_id=str(data["definition_id"]),
        definition_version=str(data["definition_version"]),
        schema_version=str(data["schema_version"]),
        default_budget=default_budget,
        maximum_budget=maximum_budget,
        allowed_actions=normalized,
        no_progress_turn_limit=int(data["no_progress_turn_limit"]),
        default_rollout_state=str(data["default_rollout_state"]),
        rollout_states=tuple(data["rollout_states"]),
    )
    policy.validate_budget(default_budget)
    if policy.no_progress_turn_limit < 1:
        raise ValueError("no_progress_turn_limit must be at least one")
    if policy.rollout_states != (
        "off",
        "shadow",
        "advisory",
        "canonical_advisory",
    ):
        raise ValueError("assistance rollout states are not the supported sequence")
    if policy.default_rollout_state != "off":
        raise ValueError("experimental assistance must default to off")
    return policy
