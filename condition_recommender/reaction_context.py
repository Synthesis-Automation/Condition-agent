"""Graph-based relations for advisory forward/retro condition context.

Generated reactions are hypotheses. Their relation to a supplied reaction never
changes that reaction's recommendation eligibility or evidence counts.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Literal

from reactive_taxonomy import canonical_molecule_collection, featurize_reaction
from reactive_taxonomy.shared_reaction_core import (
    SharedReactionCore,
    build_shared_reaction_core,
    compare_reaction_cores,
)


def load_reaction_context_policy() -> dict[str, Any]:
    """Load validated, versioned interactive search limits."""
    policy = json.loads(
        (
            Path(__file__).with_name("definitions") / "reaction_context.v1.json"
        ).read_text(encoding="utf-8")
    )
    limits = {
        "forward_operators",
        "forward_assignments",
        "forward_outcomes",
        "forward_products",
        "retro_templates_per_tier",
        "retro_validations_per_tier",
        "retro_strategies",
        "retro_realizations",
        "alternatives",
        "recipes",
    }
    if (
        set(policy) != limits | {"definition_id", "schema_version"}
        or policy["definition_id"] != "reaction_context.v1"
        or policy["schema_version"] != "1.0"
    ):
        raise ValueError("invalid reaction context definition")
    if any(type(policy[k]) is not int or policy[k] < 1 for k in limits):
        raise ValueError("reaction context limits must be positive integers")
    return policy


def reaction_sides(reaction: str) -> tuple[str, str]:
    """Validate and canonicalize complete sides without inventing missing inputs."""
    parts = reaction.split(">")
    if len(parts) != 3 or not parts[0] or not parts[2]:
        raise ValueError("reaction context requires complete reactants and products")
    sides = tuple(canonical_molecule_collection(p) for p in (parts[0], parts[2]))
    if not all(sides) or (parts[1] and not canonical_molecule_collection(parts[1])):
        raise ValueError("reaction context contains invalid molecular structures")
    return sides[0], sides[1]


def context_core(reaction: str) -> SharedReactionCore:
    """Keep unresolved/conflicting correspondence unavailable for comparison."""
    analysis = featurize_reaction(reaction)
    return build_shared_reaction_core(
        reaction,
        asdict(analysis.reaction_signature) if analysis.reaction_signature else {},
        asdict(analysis.reaction_core) if analysis.reaction_core else {},
    )


@dataclass(frozen=True)
class PlanningRelation:
    """Relation of a generated proposal to the original requested graph edits."""

    kind: Literal[
        "supplied_reaction",
        "precursor_alternative",
        "route_alternative",
        "different_product",
        "unresolved",
    ]
    inputs_changed: bool
    core_level: str | None
    differences: tuple[str, ...]
    reasons: tuple[str, ...]
    advisory_only: bool = True

    def to_dict(self) -> dict[str, Any]:
        """Serialize the relation and its uncertainty."""
        return asdict(self)


def classify_planned_reaction(
    query: str,
    proposal: str,
    *,
    query_core: SharedReactionCore | None = None,
) -> PlanningRelation:
    """Require complete graph evidence before calling a proposal a core analogue."""
    inputs, product = reaction_sides(query)
    proposed_inputs, proposed_product = reaction_sides(proposal)
    changed = inputs != proposed_inputs
    if product != proposed_product:
        return PlanningRelation(
            "different_product", changed, None, (), ("PRODUCT_GRAPH_DIFFERS",)
        )
    first = query_core if query_core is not None else context_core(query)
    second = context_core(proposal)
    if not first.levels or not second.levels:
        return PlanningRelation(
            "unresolved",
            changed,
            None,
            (),
            tuple(
                dict.fromkeys(
                    (
                        "CORE_COMPARISON_UNAVAILABLE",
                        *first.unavailable_reasons,
                        *second.unavailable_reasons,
                        *first.warnings,
                        *second.warnings,
                    )
                )
            ),
        )
    comparison = compare_reaction_cores(first, second)
    kind = (
        ("precursor_alternative" if changed else "supplied_reaction")
        if comparison.eligible
        else "route_alternative"
    )
    return PlanningRelation(
        kind, changed, comparison.level, comparison.differences, comparison.reasons
    )
