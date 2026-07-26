"""Deterministic normalization of publication and patent references."""

from __future__ import annotations

import hashlib
import re
import unicodedata
from typing import Optional

from ..models import ReferenceIdentity

_DOI_RE = re.compile(r"\b10\.\d{4,9}/[-._;()/:a-z0-9]+\b", re.IGNORECASE)
_PATENT_RE = re.compile(
    r"\b(WO|US|EP|CN|JP|KR)\s*[-/]?\s*(\d{5,})(?:\s*[-/]?\s*([A-Z]\d?))?\b",
    re.IGNORECASE,
)
_YEAR_RE = re.compile(r"\b(?:18|19|20|21)\d{2}\b")
_TRAILING_DOI_PUNCTUATION = ".,;:)]}"


def _normalized_text(value: str) -> str:
    text = unicodedata.normalize("NFKC", value)
    text = (
        text.replace("\u2010", "-")
        .replace("\u2011", "-")
        .replace("\u2012", "-")
        .replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2212", "-")
    )
    return " ".join(text.casefold().split()).strip()


def _digest(identity: str) -> str:
    return "REF1:" + hashlib.sha256(identity.encode("utf-8")).hexdigest()


def _extract_doi(normalized: str) -> Optional[str]:
    match = _DOI_RE.search(normalized)
    if match is None:
        return None
    return match.group(0).rstrip(_TRAILING_DOI_PUNCTUATION).casefold()


def _extract_patent_number(normalized: str) -> Optional[str]:
    match = _PATENT_RE.search(normalized)
    if match is None:
        return None
    country, number, kind = match.groups()
    return f"{country.upper()}{number}{(kind or '').upper()}"


def normalize_reference(reference: str) -> ReferenceIdentity:
    """Return a stable identity without treating citation text as chemistry."""
    raw = str(reference or "").strip()
    normalized = _normalized_text(raw)
    if not normalized:
        return ReferenceIdentity(
            reference_id="",
            raw_reference=raw,
            normalized_citation="",
            resolution_status="missing",
            warnings=("MISSING_REFERENCE",),
        )

    doi = _extract_doi(normalized)
    patent_number = _extract_patent_number(normalized)
    year_match = _YEAR_RE.search(normalized)
    publication_year = int(year_match.group(0)) if year_match else None
    if doi:
        identity = f"doi:{doi}"
        status = "doi"
    elif patent_number:
        identity = f"patent:{patent_number}"
        status = "patent_number"
    else:
        identity = f"citation:{normalized}"
        status = "bibliographic_text"
    warnings = (
        ("REFERENCE_ID_FROM_NORMALIZED_TEXT",)
        if status == "bibliographic_text"
        else ()
    )
    return ReferenceIdentity(
        reference_id=_digest(identity),
        raw_reference=raw,
        normalized_citation=normalized,
        doi=doi,
        patent_number=patent_number,
        publication_year=publication_year,
        resolution_status=status,
        warnings=warnings,
    )


__all__ = ["normalize_reference"]
