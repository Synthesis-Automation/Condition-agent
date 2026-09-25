"""Direct supplied-recipe compatibility assessment regressions."""

from condition_recommender import assess_reaction_recipe


def test_assess_reaction_recipe_uses_structural_signature() -> None:
    assessment = assess_reaction_recipe(
        "CCBr.N>>CCN",
        {"temperature_c": 25.0, "atmosphere": "nitrogen"},
    )

    assert assessment.compatible
    assert assessment.definition_id == "compatibility.v1"


def test_assess_reaction_recipe_retains_unresolved_without_inventing_conflict() -> None:
    assessment = assess_reaction_recipe(
        "CC>>CC",
        {"temperature_c": 25.0},
    )

    assert not assessment.compatible
    assert assessment.hard_conflicts == ()
    assert assessment.status == "unknown"
    assert assessment.analysis_status == "unsupported_or_unresolved"
    assert assessment.unresolved_requirements == ("VERIFIED_REACTION_SIGNATURE_REQUIRED",)
    assert assessment.schema_version == "reaction_recipe_assessment.v2"


def test_invalid_reaction_is_not_a_recipe_conflict() -> None:
    assessment = assess_reaction_recipe("not_a_reaction", {})
    assert assessment.status == "invalid_input"
    assert not assessment.hard_conflicts
    assert not assessment.compatible


def test_conflicting_mapping_keeps_warnings_and_does_not_force_recipe_conflict() -> None:
    assessment = assess_reaction_recipe("[CH3:1][Br:2].[NH3:1]>>[CH3:1][NH2:1]", {})
    assert assessment.status in {"unknown", "invalid_input"}
    assert not assessment.hard_conflicts
    assert assessment.analysis_warnings


def test_actual_structural_recipe_conflict_is_still_rejected() -> None:
    assessment = assess_reaction_recipe(
        "BrCCc1ccc(C=O)cc1.N>>NCCc1ccc(C=O)cc1",
        {"oxidants": [{"identity_status": "resolved", "role_status": "assigned", "primary_role": "oxidant", "roles": [{"role_id": "oxidant"}]}]},
    )
    assert assessment.status == "conflict"
    assert assessment.hard_conflicts
    assert not assessment.compatible
