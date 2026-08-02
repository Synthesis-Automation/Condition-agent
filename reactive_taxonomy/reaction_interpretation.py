"""Optional reaction-pattern interpretation over structural observations.

This layer consumes normalized graph evidence. It cannot propose atom
correspondence, edits, products, or partner roles.
"""

from __future__ import annotations

from .reaction_models import (
    ReactionInterpretation,
    ReactionObservation,
    ReactionPatternMatch,
)
from .reaction_patterns import match_reaction_patterns, select_primary_pattern


def _co_occurring_patterns(
    pattern_matches: tuple[ReactionPatternMatch, ...],
    observation: ReactionObservation,
) -> tuple[str, ...]:
    """Select non-overlapping synthesis-pattern interpretations."""
    selected = []
    occupied: set[int] = set()
    synthesis_matches = tuple(
        sorted(
            (
                match
                for match in pattern_matches
                if match.tier == "synthesis" and match.display_label
            ),
            key=lambda match: (
                not any(
                    observation.edits[index].edit_type == "formed"
                    for index in match.matched_edit_indices
                ),
                -match.specificity,
                -match.display_importance,
                match.pattern_id,
            ),
        )
    )
    for match in synthesis_matches:
        indices = set(match.matched_edit_indices)
        if (
            not indices
            or occupied.intersection(indices)
        ):
            continue
        selected.append(match.pattern_id)
        occupied.update(indices)
    return tuple(selected) if len(selected) > 1 else ()


def interpret_reaction(
    observation: ReactionObservation,
    *,
    label_style: str = "unicode",
) -> ReactionInterpretation:
    """Add optional pattern and family annotations to an observation."""
    del label_style
    pattern_matches = match_reaction_patterns(observation)
    co_occurring_pattern_ids = _co_occurring_patterns(
        pattern_matches,
        observation,
    )
    compatible_families = tuple(
        sorted(
            {
                family
                for match in pattern_matches
                for family in match.compatible_named_families
            }
        )
    )
    strong_synthesis_families = tuple(
        sorted(
            {
                family
                for match in pattern_matches
                if match.tier == "synthesis" and match.confidence >= 0.8
                and not match.requires_condition_evidence
                for family in match.compatible_named_families
            }
        )
    )
    return ReactionInterpretation(
        pattern_matches=pattern_matches,
        primary_pattern_id=(
            co_occurring_pattern_ids[0]
            if co_occurring_pattern_ids
            else select_primary_pattern(pattern_matches)
        ),
        co_occurring_pattern_ids=co_occurring_pattern_ids,
        compatible_named_families=compatible_families,
        named_family=(
            strong_synthesis_families[0]
            if len(strong_synthesis_families) == 1
            else None
        ),
        evidence_quality=(
            observation.evidence_quality if pattern_matches else "unresolved"
        ),
    )


__all__ = ["interpret_reaction"]
