"""Declarative, structure-grounded condition hypotheses for forward ranking."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Tuple


_PATH = Path(__file__).with_name("definitions") / "condition_profiles.v1.json"
CONDITION_PROFILE_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ForwardConditionProfile:
    """Coarse user-supplied reaction environment, not a resolved recipe."""

    strategy: str = "unspecified"
    redox_mode: str = "neutral"
    medium: str = "neutral"
    catalyst_family: str = "unspecified"
    source: str = "user_selected"
    schema_version: str = CONDITION_PROFILE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardConditionProfileEvidence:
    """Transparent effect of a coarse condition profile on one candidate."""

    evaluated: bool
    profile: ForwardConditionProfile
    score_adjustment: float
    matched_rules: Tuple[str, ...] = ()
    cautions: Tuple[str, ...] = ()
    definition_id: str = "forward_condition_profiles.v1"
    definition_version: str = "1.0"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@lru_cache(maxsize=1)
def load_condition_profile_definition() -> dict[str, Any]:
    """Load and validate the public profile catalog and score adjustments."""

    value = json.loads(_PATH.read_text(encoding="utf-8"))
    if value.get("definition_id") != "forward_condition_profiles.v1":
        raise ValueError("unexpected forward condition-profile definition")
    if value.get("schema_version") != CONDITION_PROFILE_SCHEMA_VERSION:
        raise ValueError("unsupported forward condition-profile schema")
    for field in ("strategies", "redox_modes", "media", "catalyst_families"):
        entries = value.get(field)
        if not isinstance(entries, list) or not entries:
            raise ValueError(f"condition-profile {field} must be nonempty")
        ids = [str(item.get("id") or "") for item in entries]
        if any(not item for item in ids) or len(ids) != len(set(ids)):
            raise ValueError(f"condition-profile {field} IDs must be unique")
    adjustments = value.get("score_adjustments")
    if not isinstance(adjustments, dict) or any(
        not isinstance(amount, (int, float)) for amount in adjustments.values()
    ):
        raise ValueError("condition-profile score adjustments are invalid")
    cap = float(value.get("maximum_absolute_adjustment") or 0.0)
    if not 0.0 < cap <= 1.0:
        raise ValueError("condition-profile adjustment cap must be in (0, 1]")
    return value


def condition_profile_catalog() -> dict[str, Any]:
    """Return the versioned UI-safe profile catalog."""

    value = load_condition_profile_definition()
    return {
        key: value[key]
        for key in (
            "definition_id",
            "schema_version",
            "strategies",
            "redox_modes",
            "media",
            "catalyst_families",
        )
    }


def normalize_condition_profile(
    value: Mapping[str, Any] | ForwardConditionProfile | None,
) -> ForwardConditionProfile:
    """Validate and normalize a profile against the declarative catalog."""

    if value is None:
        return ForwardConditionProfile()
    if isinstance(value, ForwardConditionProfile):
        profile = value
    elif isinstance(value, Mapping):
        profile = ForwardConditionProfile(
            strategy=str(value.get("strategy") or "unspecified"),
            redox_mode=str(value.get("redox_mode") or "neutral"),
            medium=str(value.get("medium") or "neutral"),
            catalyst_family=str(value.get("catalyst_family") or "unspecified"),
        )
    else:
        raise TypeError("condition profile must be a mapping")
    definition = load_condition_profile_definition()
    allowed = {
        "strategy": {item["id"] for item in definition["strategies"]},
        "redox_mode": {item["id"] for item in definition["redox_modes"]},
        "medium": {item["id"] for item in definition["media"]},
        "catalyst_family": {
            item["id"] for item in definition["catalyst_families"]
        },
    }
    for field, values in allowed.items():
        if getattr(profile, field) not in values:
            raise ValueError(f"unsupported condition-profile {field}")
    if (
        profile.strategy != "transition_metal_catalysis"
        and profile.catalyst_family != "unspecified"
    ):
        raise ValueError(
            "a catalyst family requires transition_metal_catalysis strategy"
        )
    return profile


def _formed_heavy_bond(tokens: tuple[str, ...]) -> bool:
    return any(token.startswith("formed:") for token in tokens)


def _classic_coupling(tokens: tuple[str, ...]) -> bool:
    leaving_group_loss = any(
        token.startswith(("broken:Br-C:", "broken:C-Cl:", "broken:C-I:",
                          "broken:C-Br:", "broken:Cl-C:", "broken:I-C:"))
        for token in tokens
    )
    coupling_bond = any(
        token.startswith(("formed:C-C:", "formed:C-N:", "formed:C-O:",
                          "formed:C-S:"))
        for token in tokens
    )
    return leaving_group_loss and coupling_bond


def _hydrogen_gain(tokens: tuple[str, ...]) -> bool:
    return any(
        token.startswith("hydrogen_change:") and token.endswith("NONE>SINGLE")
        for token in tokens
    )


def _hydrogen_loss(tokens: tuple[str, ...]) -> bool:
    return any(
        token.startswith("hydrogen_change:") and token.endswith("SINGLE>NONE")
        for token in tokens
    )


def assess_condition_profile(
    observed_edit_tokens: tuple[str, ...],
    profile: Mapping[str, Any] | ForwardConditionProfile | None,
) -> ForwardConditionProfileEvidence:
    """Score only structurally defensible effects of a coarse condition prior."""

    normalized = normalize_condition_profile(profile)
    evaluated = normalized != ForwardConditionProfile()
    if not evaluated:
        return ForwardConditionProfileEvidence(False, normalized, 0.0)
    definition = load_condition_profile_definition()
    amounts = definition["score_adjustments"]
    rules = []
    cautions = []
    adjustment = 0.0
    classic_coupling = _classic_coupling(observed_edit_tokens)
    heavy_formation = _formed_heavy_bond(observed_edit_tokens)
    if normalized.strategy == "transition_metal_catalysis":
        if classic_coupling:
            adjustment += amounts["transition_metal_cross_coupling"]
            rules.append("transition_metal_cross_coupling")
        elif heavy_formation:
            adjustment += amounts["transition_metal_heavy_bond_formation"]
            rules.append("transition_metal_heavy_bond_formation")
        else:
            adjustment += amounts["transition_metal_no_heavy_bond_formation"]
            rules.append("transition_metal_no_heavy_bond_formation")
        if normalized.catalyst_family != "unspecified":
            cautions.append(
                "CATALYST_FAMILY_RECORDED_NOT_USED_AS_PRODUCT_PROOF"
            )
    elif normalized.strategy == "metal_free_polar" and classic_coupling:
        adjustment += amounts["metal_free_classic_coupling_caution"]
        rules.append("metal_free_classic_coupling_caution")
        cautions.append("CLASSIC_COUPLING_PATTERN_MAY_STILL_PROCEED_METAL_FREE")
    elif normalized.strategy in {"radical", "photochemical"}:
        cautions.append("MECHANISM_NOT_RESOLVED_FROM_NET_BOND_EDITS")
    if normalized.redox_mode == "oxidative":
        key = (
            "oxidative_pattern_match" if _hydrogen_loss(observed_edit_tokens)
            else "oxidative_pattern_conflict" if _hydrogen_gain(observed_edit_tokens)
            else None
        )
        if key:
            adjustment += amounts[key]
            rules.append(key)
    elif normalized.redox_mode == "reductive":
        key = (
            "reductive_pattern_match" if _hydrogen_gain(observed_edit_tokens)
            else "reductive_pattern_conflict" if _hydrogen_loss(observed_edit_tokens)
            else None
        )
        if key:
            adjustment += amounts[key]
            rules.append(key)
    if normalized.medium != "neutral":
        cautions.append("MEDIUM_RECORDED_WITHOUT_STRUCTURE_ONLY_SELECTIVITY_RULE")
    cap = float(definition["maximum_absolute_adjustment"])
    adjustment = round(max(-cap, min(cap, adjustment)), 8)
    return ForwardConditionProfileEvidence(
        evaluated=True,
        profile=normalized,
        score_adjustment=adjustment,
        matched_rules=tuple(rules),
        cautions=tuple(cautions),
    )


__all__ = [
    "CONDITION_PROFILE_SCHEMA_VERSION",
    "ForwardConditionProfile",
    "ForwardConditionProfileEvidence",
    "assess_condition_profile",
    "condition_profile_catalog",
    "load_condition_profile_definition",
    "normalize_condition_profile",
]
