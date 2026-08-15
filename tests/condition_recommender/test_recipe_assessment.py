"""Direct supplied-recipe compatibility assessment regressions."""

from condition_recommender import assess_reaction_recipe


def test_assess_reaction_recipe_uses_structural_signature() -> None:
    assessment = assess_reaction_recipe(
        "CCBr.N>>CCN",
        {"temperature_c": 25.0, "atmosphere": "nitrogen"},
    )

    assert assessment.compatible
    assert assessment.definition_id == "compatibility.v1"


def test_assess_reaction_recipe_rejects_unresolved_reaction() -> None:
    assessment = assess_reaction_recipe(
        "CC>>CC",
        {"temperature_c": 25.0},
    )

    assert not assessment.compatible
    assert assessment.hard_conflicts == ("UNRESOLVED_REACTION_FOR_RECIPE_ASSESSMENT",)
