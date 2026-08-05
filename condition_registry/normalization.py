"""Deterministic identifier normalization without legacy dependencies."""

from __future__ import annotations

import re
from typing import Optional

from .models import CONDITION_NAME_IDENTIFIER_TYPES

_CAS_RE = re.compile(r"^(\d{2,7})-(\d{2})-(\d)$")

_IDENTIFIER_NORMALIZATION_PROFILES = {
    "canonical_name": "chemical_name_v1",
    "common_name": "chemical_name_v1",
    "systematic_name": "chemical_name_v1",
    "abbreviation": "abbreviation_v1",
    "trade_name": "trade_name_v1",
    "legacy_name": "chemical_name_v1",
    "cas": "cas_v1",
    "inchi_key": "inchi_key_v1",
    "database_id": "database_id_v1",
}


def normalize_name(value: str) -> str:
    """Return the legacy broad name normalization used by v1 definitions."""
    text = str(value or "").lower().replace("·", ".")
    text = re.sub(r"\b(anhydrous|dry|aqueous|aq|solution|soln)\b", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def normalize_chemical_name(value: str) -> str:
    """Normalize a chemical name without erasing formulation distinctions."""
    text = str(value or "").casefold().replace("·", ".")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def normalize_cas(value: str) -> Optional[str]:
    text = str(value or "").strip()
    match = _CAS_RE.fullmatch(text)
    if match is None:
        return None
    digits = match.group(1) + match.group(2)
    checksum = sum(position * int(digit) for position, digit in enumerate(reversed(digits), start=1))
    return text if checksum % 10 == int(match.group(3)) else None


def normalize_identifier(value: str, identifier_type: str) -> Optional[str]:
    """Normalize one typed condition identifier for deterministic indexing."""
    if identifier_type == "cas":
        return normalize_cas(value)
    if identifier_type in CONDITION_NAME_IDENTIFIER_TYPES:
        return normalize_chemical_name(value)
    if identifier_type == "inchi_key":
        text = re.sub(r"\s+", "", str(value or "")).upper()
        return text or None
    if identifier_type == "database_id":
        text = " ".join(str(value or "").casefold().split())
        return text or None
    return None


def identifier_normalization_profile(identifier_type: str) -> Optional[str]:
    """Return the versioned normalization profile for an identifier type."""
    return _IDENTIFIER_NORMALIZATION_PROFILES.get(identifier_type)


__all__ = [
    "normalize_cas",
    "normalize_chemical_name",
    "normalize_identifier",
    "normalize_name",
    "identifier_normalization_profile",
]
