"""Shadow parity comparison between legacy and connectivity rewrites."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence, Tuple

from .connectivity_rewrite import (
    apply_connectivity_rewrite,
    connectivity_rewrite_for_grammar,
)
from .reaction_models import (
    OperatorOutcome,
    ReactionComponent,
    ReactionSiteReference,
)
from .reaction_operators import enumerate_legacy_operator_outcomes


@dataclass(frozen=True)
class ConnectivityRewriteParity:
    """Exact shadow comparison for one grammar assignment."""

    grammar_id: str
    migrated: bool
    parity: bool
    product_parity: bool
    edit_parity: bool
    warning_parity: bool
    ambiguity_parity: bool
    ordering_parity: bool
    legacy_outcomes: Tuple[OperatorOutcome, ...]
    rewrite_outcomes: Tuple[OperatorOutcome, ...]


def _products(outcomes: Sequence[OperatorOutcome]) -> tuple[tuple[str, str | None], ...]:
    return tuple(
        (outcome.outcome_id, outcome.predicted_product_smiles)
        for outcome in outcomes
    )


def _edits(outcomes: Sequence[OperatorOutcome]) -> tuple[tuple[str, tuple], ...]:
    return tuple(
        (outcome.outcome_id, outcome.predicted_bond_changes)
        for outcome in outcomes
    )


def _warnings(outcomes: Sequence[OperatorOutcome]) -> tuple[tuple[str, tuple], ...]:
    return tuple((outcome.outcome_id, outcome.warnings) for outcome in outcomes)


def compare_connectivity_rewrite(
    grammar: Mapping[str, Any],
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
) -> ConnectivityRewriteParity:
    """Dual-run one assignment while leaving legacy results authoritative."""
    legacy = enumerate_legacy_operator_outcomes(
        dict(grammar), dict(assignment), tuple(components)
    )
    rewrite = apply_connectivity_rewrite(grammar, assignment, components)
    migrated = (
        connectivity_rewrite_for_grammar(str(grammar.get("id") or "")) is not None
    )
    product_parity = migrated and _products(legacy) == _products(rewrite)
    edit_parity = migrated and _edits(legacy) == _edits(rewrite)
    warning_parity = migrated and _warnings(legacy) == _warnings(rewrite)
    ambiguity_parity = migrated and len(legacy) == len(rewrite)
    ordering_parity = migrated and tuple(
        outcome.outcome_id for outcome in legacy
    ) == tuple(outcome.outcome_id for outcome in rewrite)
    parity = (
        product_parity
        and edit_parity
        and warning_parity
        and ambiguity_parity
        and ordering_parity
    )
    return ConnectivityRewriteParity(
        grammar_id=str(grammar.get("id") or ""),
        migrated=migrated,
        parity=parity,
        product_parity=product_parity,
        edit_parity=edit_parity,
        warning_parity=warning_parity,
        ambiguity_parity=ambiguity_parity,
        ordering_parity=ordering_parity,
        legacy_outcomes=legacy,
        rewrite_outcomes=rewrite,
    )


__all__ = ["ConnectivityRewriteParity", "compare_connectivity_rewrite"]
