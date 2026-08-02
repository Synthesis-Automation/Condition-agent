"""Stable review serialization for optional reaction-pattern interpretation."""

from __future__ import annotations

from typing import Any, Dict, Mapping


REACTION_PATTERN_REVIEW_FIELDS = (
    "primary_reaction_pattern",
    "primary_reaction_pattern_count",
    "reaction_pattern_matches",
    "co_occurring_reaction_patterns",
    "identified_reaction_type",
    "compatible_reaction_types",
    "reaction_pattern_confidence",
    "reaction_pattern_requires_condition_evidence",
)


def flattened_reaction_pattern_fields(
    interpretation: Mapping[str, Any] | None,
    *,
    named_family: Any = None,
    separator: str = "|",
) -> Dict[str, str]:
    """Flatten canonical pattern interpretation without recomputing chemistry."""
    value = interpretation if isinstance(interpretation, Mapping) else {}
    matches = tuple(
        match
        for match in value.get("pattern_matches") or ()
        if isinstance(match, Mapping) and match.get("pattern_id")
    )
    primary_id = str(value.get("primary_pattern_id") or "")
    primary = next(
        (match for match in matches if str(match.get("pattern_id")) == primary_id),
        None,
    )
    compatible = tuple(
        str(family)
        for family in value.get("compatible_named_families") or ()
        if family
    )
    return {
        "primary_reaction_pattern": primary_id,
        "primary_reaction_pattern_count": (
            str(primary.get("occurrence_count") or 1)
            if primary is not None
            else ""
        ),
        "reaction_pattern_matches": separator.join(
            str(match["pattern_id"]) for match in matches
        ),
        "co_occurring_reaction_patterns": separator.join(
            str(pattern_id)
            for pattern_id in value.get("co_occurring_pattern_ids") or ()
            if pattern_id
        ),
        "identified_reaction_type": str(
            named_family or value.get("named_family") or ""
        ),
        "compatible_reaction_types": separator.join(compatible),
        "reaction_pattern_confidence": (
            str(primary.get("confidence")) if primary is not None else ""
        ),
        "reaction_pattern_requires_condition_evidence": (
            str(bool(primary.get("requires_condition_evidence")))
            if primary is not None
            else ""
        ),
    }


__all__ = [
    "REACTION_PATTERN_REVIEW_FIELDS",
    "flattened_reaction_pattern_fields",
]
