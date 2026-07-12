"""Deterministic identifier normalization without legacy dependencies."""

from __future__ import annotations

import re
from typing import Optional

_CAS_RE = re.compile(r"^(\d{2,7})-(\d{2})-(\d)$")


def normalize_name(value: str) -> str:
    text = str(value or "").lower().replace("·", ".")
    text = re.sub(r"\b(anhydrous|dry|aqueous|aq|solution|soln)\b", " ", text)
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


__all__ = ["normalize_cas", "normalize_name"]
