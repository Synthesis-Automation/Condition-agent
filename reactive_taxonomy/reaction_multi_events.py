"""Exact reconstruction of balanced, unmapped multi-event reactions."""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, Iterable, Sequence, Tuple

from .reaction_models import ReactionComponent, ReactionSiteReference
from .reaction_operators import apply_operator_sequence


RawCandidate = Tuple[Dict[str, Any], Dict[str, ReactionSiteReference]]


def _operation_key(candidate: RawCandidate) -> Tuple[Any, ...]:
    grammar, assignment = candidate
    operator = grammar["operator"]
    electrophile = assignment[str(operator["electrophile_role"])]
    partner = assignment[str(operator["partner_role"])]
    anchor_role = "anchor" if "anchor" in electrophile.atom_roles else "center"
    return (
        grammar["id"],
        electrophile.component_index,
        int(electrophile.atom_roles[anchor_role][0]),
        partner.component_index,
        int(partner.atom_roles["center"][0]),
    )


def _interpretation_key(candidates: Sequence[RawCandidate]) -> Tuple[Any, ...]:
    return tuple(
        sorted(
            (
                grammar["id"],
                tuple(
                    sorted(
                        site.canonical_signature for site in assignment.values()
                    )
                ),
            )
            for grammar, assignment in candidates
        )
    )


def exact_multi_event_reconstructions(
    raw_candidates: Iterable[RawCandidate],
    reactants: Tuple[ReactionComponent, ...],
    observed_products: set[str],
    *,
    max_events: int = 4,
    max_combinations: int = 5000,
) -> Tuple[Tuple[RawCandidate, ...], ...]:
    """Return composite operator sets that exactly reconstruct one product."""
    eligible = tuple(
        sorted(
            (
                candidate
                for candidate in raw_candidates
                if candidate[0].get("operator", {}).get("id")
                == "center_replacement"
            ),
            key=_operation_key,
        )
    )
    if len(eligible) < 2 or not observed_products:
        return ()
    exact = []
    attempted = 0
    upper = min(max_events, len(eligible))
    for event_count in range(2, upper + 1):
        for selected in combinations(eligible, event_count):
            attempted += 1
            if attempted > max_combinations:
                return tuple(exact)
            predicted = apply_operator_sequence(selected, reactants)
            if predicted in observed_products:
                exact.append(tuple(selected))
        if exact:
            break
    if not exact:
        return ()
    return tuple(
        sorted(
            exact,
            key=lambda selected: (
                _interpretation_key(selected),
                tuple(_operation_key(candidate) for candidate in selected),
            ),
        )
    )


def equivalent_multi_event_interpretations(
    reconstructions: Sequence[Sequence[RawCandidate]],
) -> bool:
    """Return whether exact alternatives differ only by equivalent sites."""
    return len({_interpretation_key(selected) for selected in reconstructions}) <= 1


__all__ = [
    "RawCandidate",
    "equivalent_multi_event_interpretations",
    "exact_multi_event_reconstructions",
]
