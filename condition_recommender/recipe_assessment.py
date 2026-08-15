"""Assess one supplied resolved recipe against a fully specified reaction."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, Mapping

from reactive_taxonomy import featurize_reaction

from .compatibility import CompatibilityAssessment, assess_recipe_compatibility


def assess_reaction_recipe(
    reaction_smiles: str,
    recipe: Mapping[str, Any],
) -> CompatibilityAssessment:
    """Apply canonical compatibility rules to a proposed reaction and recipe.

    ``recipe`` is expected to be a resolved condition-registry recipe. This
    assesses compatibility only; it does not predict outcome or yield.
    """

    if not isinstance(recipe, Mapping):
        raise TypeError("recipe must be a resolved recipe mapping")
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid or analysis.reaction_signature is None:
        return CompatibilityAssessment(
            compatible=False,
            score=0.0,
            hard_conflicts=("UNRESOLVED_REACTION_FOR_RECIPE_ASSESSMENT",),
            evidence=("The proposed reaction has no verified structural signature.",),
            definition_id="compatibility.v1",
            definition_version="1.2",
        )
    return assess_recipe_compatibility(
        asdict(analysis.reaction_signature),
        recipe,
    )


__all__ = ["assess_reaction_recipe"]
