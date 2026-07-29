"""Stable flat projections from canonical typed reactivity profiles."""

from __future__ import annotations

from typing import Any, Mapping, Tuple


def profile_classes(partner: Mapping[str, Any]) -> Tuple[str, str]:
    """Return accessibility and activation classes for review-only flat exports."""
    profile = partner.get("reactivity_profile") or {}
    if not isinstance(profile, Mapping):
        return "", ""
    steric = profile.get("steric") or {}
    electronic = profile.get("electronic") or {}
    steric_value = (
        str(steric.get("accessibility_class") or "")
        if isinstance(steric, Mapping)
        else ""
    )
    electronic_value = (
        str(electronic.get("activation_class") or "")
        if isinstance(electronic, Mapping)
        else ""
    )
    return steric_value, electronic_value


__all__ = ["profile_classes"]
