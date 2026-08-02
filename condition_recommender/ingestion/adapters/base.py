"""Shared protocol and helpers for deterministic source adapters."""

from __future__ import annotations

import csv
import hashlib
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, Optional, Protocol, Sequence, Tuple

from ..models import CanonicalSourceObservation, SourceProvenance


_MAPPED_ATOM = re.compile(r":\d+\]")


class SourceAdapter(Protocol):
    """Executable adapter selected from the explicit in-process registry."""

    adapter_id: str
    adapter_version: str
    corpus_id: str
    required_columns: Tuple[str, ...]

    def iter_observations(
        self, path: Path, *, source_sha256: str
    ) -> Iterator[CanonicalSourceObservation]: ...


def read_headers(path: Path) -> Tuple[str, ...]:
    """Read a CSV header without loading source rows."""
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return tuple(str(value).strip() for value in next(csv.reader(handle), ()))


def validate_headers(path: Path, required: Sequence[str]) -> None:
    """Reject a file whose exact source contract is incomplete."""
    headers = set(read_headers(path))
    missing = sorted(set(required).difference(headers))
    if missing:
        raise ValueError(
            f"{path.name} is missing required columns: {', '.join(missing)}"
        )


def clean_text(value: Any) -> str:
    """Return source text with surrounding whitespace removed."""
    return str(value or "").strip()


def optional_float(value: Any) -> tuple[Optional[float], Optional[str]]:
    """Parse a finite float and return a warning code on invalid text."""
    text = clean_text(value)
    if not text:
        return None, None
    try:
        result = float(text)
    except (TypeError, ValueError):
        return None, "INVALID_NUMERIC_VALUE"
    if not math.isfinite(result):
        return None, "NONFINITE_NUMERIC_VALUE"
    return result, None


def optional_int(value: Any) -> tuple[Optional[int], Optional[str]]:
    """Parse an integer-valued source field."""
    number, warning = optional_float(value)
    if warning or number is None:
        return None, warning
    if not number.is_integer():
        return None, "INVALID_INTEGER_VALUE"
    return int(number), None


def split_identifiers(value: Any) -> Tuple[str, ...]:
    """Split source fields that explicitly contain identifier lists."""
    return tuple(
        item.strip() for item in re.split(r"[,;|]", clean_text(value)) if item.strip()
    )


def supplied_mapping_status(reaction_smiles: str) -> str:
    """Describe only whether atom-map syntax was supplied, without validating it."""
    if not reaction_smiles:
        return "missing_structure"
    return (
        "supplied_unvalidated"
        if _MAPPED_ATOM.search(reaction_smiles)
        else "not_supplied"
    )


def source_provenance(
    *,
    adapter: SourceAdapter,
    path: Path,
    source_sha256: str,
    row_number: int,
    record_id: str,
    source_groups: Optional[Dict[str, str]] = None,
    reference: str = "",
) -> SourceProvenance:
    """Build consistent provenance without using machine-specific absolute paths."""
    return SourceProvenance(
        corpus_id=adapter.corpus_id,
        release_id=path.stem,
        adapter_id=adapter.adapter_id,
        adapter_version=adapter.adapter_version,
        source_file=path.name,
        source_file_sha256=source_sha256,
        source_row_number=row_number,
        source_record_id=record_id,
        source_groups=dict(source_groups or {}),
        reference=reference,
    )


def observation_id(
    *, adapter_id: str, source_sha256: str, row_number: int, record_id: str
) -> str:
    """Identify one immutable source-file row independently of its local path."""
    token = "\x1f".join(
        (adapter_id, source_sha256, str(row_number), clean_text(record_id))
    )
    return "SOBS1:" + hashlib.sha256(token.encode("utf-8")).hexdigest()


def raw_fields(row: Mapping[str, Any]) -> Dict[str, Any]:
    """Preserve every named source field for lossless audit."""
    return {str(key): value for key, value in row.items() if key is not None}


__all__ = [
    "SourceAdapter",
    "clean_text",
    "observation_id",
    "optional_float",
    "optional_int",
    "raw_fields",
    "read_headers",
    "source_provenance",
    "split_identifiers",
    "supplied_mapping_status",
    "validate_headers",
]
