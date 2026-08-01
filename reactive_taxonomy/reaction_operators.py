"""Grammar-independent public graph-operator registry.

The compatibility rewrite loader still retains grammar-to-role adapters, but
this public surface exposes only executable structural operator definitions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence, Tuple

from .connectivity_rewrite import (
    CompiledRewriteVariant,
    apply_reaction_operator as _apply_reaction_operator,
    load_connectivity_rewrites,
)
from .reaction_models import (
    ReactionComponent,
    ReactionSiteReference,
    RewriteOutcome,
)


@dataclass(frozen=True)
class ReactionOperator:
    """Executable structural transformation without grammar semantics."""

    operator_id: str
    template: str
    variants: Tuple[CompiledRewriteVariant, ...]
    schema_version: str
    instruction_set_version: str
    site_interface_schema_version: str


def load_reaction_operators() -> Tuple[ReactionOperator, ...]:
    """Return structural views of all registered connectivity operators."""
    return tuple(
        ReactionOperator(
            operator_id=rewrite.rewrite_id,
            template=rewrite.template,
            variants=rewrite.variants,
            schema_version=rewrite.schema_version,
            instruction_set_version=rewrite.instruction_set_version,
            site_interface_schema_version=rewrite.site_interface_schema_version,
        )
        for rewrite in load_connectivity_rewrites()
    )


def get_reaction_operator(operator_id: str) -> ReactionOperator | None:
    """Return a registered structural operator by ID."""
    return next(
        (
            operator
            for operator in load_reaction_operators()
            if operator.operator_id == operator_id
        ),
        None,
    )


def execute_reaction_operator(
    operator_id: str,
    assignment: Mapping[str, ReactionSiteReference],
    components: Sequence[ReactionComponent],
) -> Tuple[RewriteOutcome, ...]:
    """Execute a structural operator against generic site-interface bindings."""
    return _apply_reaction_operator(operator_id, assignment, components)


__all__ = [
    "ReactionOperator",
    "execute_reaction_operator",
    "get_reaction_operator",
    "load_reaction_operators",
]
