"""Optional reaction-pattern interpretation over structural observations.

This layer consumes normalized graph evidence. It cannot propose atom
correspondence, edits, products, or partner roles.
"""

from __future__ import annotations

from .reaction_models import ReactionInterpretation, ReactionObservation
from .reaction_patterns import match_reaction_patterns, select_primary_pattern


def interpret_reaction(
    observation: ReactionObservation,
    *,
    label_style: str = "unicode",
) -> ReactionInterpretation:
    """Add optional pattern and family annotations to an observation."""
    del label_style
    pattern_matches = match_reaction_patterns(observation)
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
                for family in match.compatible_named_families
            }
        )
    )
    return ReactionInterpretation(
        pattern_matches=pattern_matches,
        primary_pattern_id=select_primary_pattern(pattern_matches),
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
