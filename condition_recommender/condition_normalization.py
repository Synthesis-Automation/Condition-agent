"""Deterministic normalization of source condition identifiers and numbers."""

from __future__ import annotations

import re
from typing import Any, Optional, Tuple

from condition_registry.normalization import normalize_cas

from .models import ConditionIdentity

def is_valid_cas(value: str) -> bool:
    """Validate CAS syntax and check digit."""
    return normalize_cas(value) is not None


def split_identifiers(value: Any) -> Tuple[str, ...]:
    """Return stable, de-duplicated comma/semicolon-separated identifiers."""
    raw = re.split(r"[,;|]", str(value or ""))
    return tuple(sorted({item.strip() for item in raw if item.strip()}))


def normalize_cas_list(value: Any) -> Tuple[str, ...]:
    """Keep syntactically valid CAS-like identifiers in canonical order."""
    return tuple(item for item in split_identifiers(value) if is_valid_cas(item))


def optional_float(value: Any) -> Optional[float]:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    return number


def normalize_conditions(row: dict[str, Any]) -> ConditionIdentity:
    return ConditionIdentity(
        catalyst_cas=normalize_cas_list(row.get("catalyst_cas")),
        reagent_cas=normalize_cas_list(row.get("reagent_cas")),
        solvent_cas=normalize_cas_list(row.get("solvent_cas")),
    )


__all__ = ["is_valid_cas", "normalize_cas_list", "normalize_conditions", "optional_float", "split_identifiers"]
