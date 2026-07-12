"""Deterministic normalization of source condition identifiers and numbers."""

from __future__ import annotations

import re
from typing import Any, Optional, Tuple

from .models import ConditionIdentity

_CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


def is_valid_cas(value: str) -> bool:
    """Validate CAS syntax and check digit."""
    match = _CAS_RE.fullmatch(str(value or "").strip())
    if match is None:
        return False
    body, check = value.rsplit("-", 1)
    digits = body.replace("-", "")
    checksum = sum(position * int(digit) for position, digit in enumerate(reversed(digits), start=1))
    return checksum % 10 == int(check)


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
