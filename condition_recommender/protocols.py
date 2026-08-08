"""Compose registry-owned protocol drafts from reaction observations."""

from __future__ import annotations

from typing import Any, Mapping

from condition_registry import (
    ProtocolReactionMaterialInput,
    ResolvedConditionRecipe,
    SynthesisProtocolDraft,
    build_synthesis_protocol_draft,
)
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


def protocol_draft_for_reaction(
    recipe: ResolvedConditionRecipe | Mapping[str, Any],
    reaction_smiles: str,
) -> SynthesisProtocolDraft:
    """Add graph-parsed reaction materials to a canonical condition recipe."""
    parsed = parse_reaction_smiles(
        reaction_smiles,
        include_molecular_interpretation=False,
    )
    materials = ()
    if parsed.valid:
        materials = tuple(
            ProtocolReactionMaterialInput(
                side=component.side,  # type: ignore[arg-type]
                component_index=component.component_index,
                smiles=component.input_smiles,
                canonical_smiles=component.canonical_smiles or None,
            )
            for component in (*parsed.reactants, *parsed.agents, *parsed.products)
        )
    return build_synthesis_protocol_draft(
        recipe,
        reaction_smiles=reaction_smiles or None,
        reaction_materials=materials,
    )


__all__ = ["protocol_draft_for_reaction"]
