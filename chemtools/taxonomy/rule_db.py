"""
Helpers for integrating taxonomy reaction types with rule database assets.

This module resolves a detected (or user-provided) reaction family label to a
rule database identifier stored in taxonomy v2 metadata:

    "metadata": {"rule_db_v2": "<db_stem_or_filename>"}
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any, Iterable, List, Optional

from . import reaction_catalog
from ..featurizers.analysis.reactions import canonical_family_label


_RULE_DB_V2_METADATA_KEY = "rule_db_v2"


def _slugify(value: str) -> str:
    return re.sub(r"[^0-9a-z]+", "_", value.lower()).strip("_")


def _dedupe(items: Iterable[str]) -> List[str]:
    seen: set[str] = set()
    out: List[str] = []
    for item in items:
        if not item or item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def _reaction_type_candidates(label: str) -> List[str]:
    raw = (label or "").strip()
    if not raw:
        return []

    variants: List[str] = [raw]

    lower = raw.lower()
    if lower != raw:
        variants.append(lower)

    slug = _slugify(raw)
    if slug and slug not in variants:
        variants.append(slug)

    # Common alternates for alias keys: underscores vs spaces/hyphens.
    for token in (raw, lower, slug):
        if not token:
            continue
        if "_" in token:
            variants.append(token.replace("_", " "))
            variants.append(token.replace("_", "-"))
        if "-" in token:
            variants.append(token.replace("-", " "))
            variants.append(token.replace("-", "_"))

    return _dedupe(variants)


def resolve_reaction_type_id(label: str) -> Optional[str]:
    """Resolve a reaction family label to a canonical taxonomy v2 reaction id."""
    canonical = canonical_family_label(label)
    if canonical:
        return canonical
    for candidate in _reaction_type_candidates(label):
        resolved = reaction_catalog.resolve_reaction_type(candidate)
        if resolved:
            return resolved
    return None


def _normalize_rule_db_identifier(value: str) -> Optional[str]:
    text = (value or "").strip()
    if not text:
        return None
    if text.lower().endswith(".json"):
        return Path(text).stem
    return text


def resolve_rule_db_v2(family: str) -> Optional[str]:
    """
    Resolve a reaction family label to a `data/rule_db_v2` database stem.

    The returned value is a file *stem* (without `.json`) suitable for looking up
    `data/rule_db_v2/<stem>.json`.
    """
    reaction_type_id = resolve_reaction_type_id(family)
    if not reaction_type_id:
        return None

    reaction_type = reaction_catalog.get_reaction_type(reaction_type_id)
    if reaction_type is None:
        return None

    metadata: dict[str, Any] = dict(reaction_type.metadata or {})
    raw_value = metadata.get(_RULE_DB_V2_METADATA_KEY)
    if isinstance(raw_value, str):
        return _normalize_rule_db_identifier(raw_value)
    if isinstance(raw_value, list):
        for entry in raw_value:
            if isinstance(entry, str):
                resolved = _normalize_rule_db_identifier(entry)
                if resolved:
                    return resolved
        return None

    return None
