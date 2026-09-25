"""Assess one supplied resolved recipe against a fully specified reaction."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.reaction_models import ReactionAnalysis

from .compatibility import CompatibilityAssessment, assess_recipe_compatibility, load_compatibility_rules


@dataclass(frozen=True)
class ReactionRecipeAssessment(CompatibilityAssessment):
    """Reaction-level coverage plus recipe rules; false is not always a conflict.

    Inspect status and hard_conflicts: unknown/invalid_input with compatible=False
    means assessment cannot establish admission, not demonstrated incompatibility.
    """

    analysis_status: str = "supported"
    analysis_warnings: tuple[str, ...] = ()
    schema_version: str = "reaction_recipe_assessment.v2"


def assess_reaction_recipe(
    reaction_smiles: str,
    recipe: Mapping[str, Any],
) -> ReactionRecipeAssessment:
    """Apply canonical compatibility rules to a proposed reaction and recipe.

    ``recipe`` is expected to be a resolved condition-registry recipe. This
    assesses compatibility only; it does not predict outcome or yield.
    """

    if not isinstance(recipe, Mapping):
        raise TypeError("recipe must be a resolved recipe mapping")
    return _assess_analyzed_recipe(featurize_reaction(reaction_smiles), recipe)


def _assess_analyzed_recipe(analysis: ReactionAnalysis, recipe: Mapping[str, Any]) -> ReactionRecipeAssessment:
    """Reuse one unchanged analysis when inspecting several recipes."""
    if not analysis.valid or analysis.reaction_signature is None:
        rules = load_compatibility_rules()
        return ReactionRecipeAssessment(
            compatible=False,
            score=0.0,
            evidence=("The proposed reaction has no verified structural signature.",),
            definition_id=rules["definition_id"],
            definition_version=rules["schema_version"],
            status="unknown" if analysis.valid else "invalid_input",
            analysis_status="unsupported_or_unresolved" if analysis.valid else "invalid_input",
            analysis_warnings=tuple(analysis.warnings),
            unresolved_requirements=("VERIFIED_REACTION_SIGNATURE_REQUIRED",),
        )
    assessment = assess_recipe_compatibility(
        {**asdict(analysis.reaction_signature),
         "spectator_groups": tuple(asdict(group) for group in analysis.spectator_groups)},
        recipe,
    )
    return ReactionRecipeAssessment(**asdict(assessment), analysis_warnings=tuple(analysis.warnings))


__all__ = ["ReactionRecipeAssessment", "assess_reaction_recipe"]
