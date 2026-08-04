"""Validated chemist-selectable ranking profiles and custom priorities."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

from .models import ChemistRankingPreferences


_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "chemist_ranking_profiles.v1.json"
)
_DEFINITION_ID = "chemist_ranking_profiles.v1"
_SCHEMA_VERSION = "1.0"
RANKING_COMPONENTS = (
    "similarity",
    "partner_category",
    "functional_group_tolerance",
    "yield",
    "independent_support",
    "reaction_breadth",
    "dataset_diversity",
    "compatibility",
    "condition_certainty",
)


def normalize_ranking_weights(
    weights: Mapping[str, Any], *, label: str
) -> Dict[str, float]:
    """Validate and normalize one complete ranking-weight mapping."""
    if set(weights) != set(RANKING_COMPONENTS):
        missing = sorted(set(RANKING_COMPONENTS) - set(weights))
        extra = sorted(set(weights) - set(RANKING_COMPONENTS))
        raise ValueError(
            f"{label} weights do not match component vocabulary; "
            f"missing={missing}, extra={extra}"
        )
    normalized = {name: float(weights[name]) for name in RANKING_COMPONENTS}
    if any(value < 0.0 for value in normalized.values()):
        raise ValueError(f"{label} weights must be non-negative")
    total = sum(normalized.values())
    if total <= 0.0:
        raise ValueError(f"{label} weights must have positive total")
    return {name: value / total for name, value in normalized.items()}


@lru_cache(maxsize=1)
def load_chemist_ranking_profiles() -> Dict[str, Any]:
    """Load and validate versioned chemist-facing ranking presets."""
    with _DEFINITION_PATH.open("r", encoding="utf-8") as handle:
        payload = dict(json.load(handle))
    if payload.get("schema_version") != _SCHEMA_VERSION:
        raise ValueError("unsupported chemist ranking profile schema")
    if payload.get("definition_id") != _DEFINITION_ID:
        raise ValueError("unexpected chemist ranking profile definition ID")
    profiles = payload.get("profiles")
    if not isinstance(profiles, Mapping) or not profiles:
        raise ValueError("chemist ranking profiles cannot be empty")
    default_id = str(payload.get("default_profile_id") or "")
    if default_id not in profiles:
        raise ValueError("default chemist ranking profile is unavailable")
    validated = {}
    for profile_id, profile in profiles.items():
        if not isinstance(profile, Mapping):
            raise ValueError(f"invalid chemist ranking profile: {profile_id}")
        validated[str(profile_id)] = {
            **dict(profile),
            "weights": normalize_ranking_weights(
                profile.get("weights") or {},
                label=f"chemist ranking profile {profile_id}",
            ),
        }
    payload["profiles"] = validated
    tolerance = payload.get("functional_group_tolerance") or {}
    if float(tolerance.get("distance_decay") or 0.0) < 0.0:
        raise ValueError("functional-group tolerance decay cannot be negative")
    return payload


def available_ranking_profiles() -> Tuple[Dict[str, str], ...]:
    """Return stable profile metadata for CLI and GUI selection."""
    rules = load_chemist_ranking_profiles()
    return tuple(
        {
            "profile_id": profile_id,
            "label": str(profile.get("label") or profile_id),
            "description": str(profile.get("description") or ""),
        }
        for profile_id, profile in rules["profiles"].items()
    )


def resolve_ranking_preferences(
    preferences: ChemistRankingPreferences | None,
) -> ChemistRankingPreferences:
    """Resolve a named preset or complete custom weight mapping."""
    rules = load_chemist_ranking_profiles()
    requested = preferences or ChemistRankingPreferences()
    profile_id = requested.profile_id or str(rules["default_profile_id"])
    if requested.weights:
        weights = normalize_ranking_weights(
            requested.weights,
            label=f"chemist ranking override {profile_id}",
        )
        customized = True
    else:
        profile = rules["profiles"].get(profile_id)
        if profile is None:
            raise ValueError(f"unknown chemist ranking profile: {profile_id}")
        weights = dict(profile["weights"])
        customized = profile_id != str(rules["default_profile_id"])
    return ChemistRankingPreferences(
        profile_id=profile_id,
        weights=weights,
        definition_id=str(rules["definition_id"]),
        definition_version=str(rules["schema_version"]),
        customized=customized,
    )


def functional_group_distance_decay() -> float:
    """Return the declared graph-distance decay for tolerance evidence."""
    rules = load_chemist_ranking_profiles()
    return float(rules["functional_group_tolerance"]["distance_decay"])


__all__ = [
    "RANKING_COMPONENTS",
    "available_ranking_profiles",
    "functional_group_distance_decay",
    "load_chemist_ranking_profiles",
    "normalize_ranking_weights",
    "resolve_ranking_preferences",
]
