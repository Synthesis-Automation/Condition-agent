"""Direct supplied-recipe compatibility assessment regressions."""

from condition_recommender import assess_reaction_recipe
from condition_recommender import CompatibilityCoverage
from condition_recommender.compatibility import load_compatibility_rules


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
    assert assessment.coverage.assessment_status == "not_assessed"
    assert not assessment.coverage.evaluated_hard_conflict_rule_ids


def test_invalid_reaction_is_not_a_recipe_conflict() -> None:
    assessment = assess_reaction_recipe("not_a_reaction", {})
    assert assessment.status == "invalid_input"
    assert not assessment.hard_conflicts
    assert not assessment.compatible
    assert assessment.coverage.assessment_status == "not_assessed"


def test_conflicting_mapping_keeps_warnings_and_does_not_force_recipe_conflict() -> None:
    assessment = assess_reaction_recipe("[CH3:1][Br:2].[NH3:1]>>[CH3:1][NH2:1]", {})
    assert assessment.status in {"unknown", "invalid_input"}
    assert not assessment.hard_conflicts
    assert assessment.analysis_warnings
    assert assessment.coverage.assessment_status == "not_assessed"


def test_actual_structural_recipe_conflict_is_still_rejected() -> None:
    assessment = assess_reaction_recipe(
        "BrCCc1ccc(C=O)cc1.N>>NCCc1ccc(C=O)cc1",
        {"oxidants": [{"identity_status": "resolved", "role_status": "assigned", "primary_role": "oxidant", "roles": [{"role_id": "oxidant"}]}]},
    )
    assert assessment.status == "conflict"
    assert assessment.hard_conflicts
    assert not assessment.compatible
    assert assessment.coverage.assessment_status == "assessed"
    assert set(assessment.hard_conflicts) <= set(assessment.coverage.evaluated_hard_conflict_rule_ids)
    assert assessment.coverage.evaluated_soft_penalty_rule_ids == ()


def test_solvent_only_check_reports_coverage_without_changing_admission() -> None:
    recipe = {"solvents": [{
        "identity_status": "resolved", "substance_id": "cas:64-17-5",
        "role_status": "assigned", "primary_role": "solvent",
    }]}
    first = assess_reaction_recipe("CCBr.N>>CCN", recipe)
    second = assess_reaction_recipe("N.CCBr>>CCN", recipe)
    assert isinstance(first.coverage, CompatibilityCoverage)
    assert first.coverage == second.coverage
    assert first.compatible and first.score == 1.0
    assert first.status == "no_known_conflict"
    assert first.checked_requirements == ()
    assert first.coverage.schema_version == "compatibility_coverage.v1"
    assert first.coverage.capability_status == "not_covered"
    assert first.coverage.condition_identity_status == "resolved"
    rules = load_compatibility_rules()
    assert first.coverage.evaluated_hard_conflict_rule_ids == tuple(r["id"] for r in rules["hard_conflicts"])
    assert first.coverage.evaluated_soft_penalty_rule_ids == tuple(r["id"] for r in rules["soft_penalties"])
    assert first.coverage.score_meaning == rules["coverage_policy"]["score_meaning"]


def test_capability_coverage_distinguishes_supplied_and_missing_evidence() -> None:
    reaction = "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]"
    supported = assess_reaction_recipe(reaction, {"atmosphere": "hydrogen"})
    missing = assess_reaction_recipe(reaction, {"atmosphere": "nitrogen"})
    assert supported.coverage.capability_status == "supported"
    assert supported.checked_requirements
    assert missing.coverage.capability_status == "unresolved"
    assert missing.unresolved_requirements
    assert missing.compatible and not missing.hard_conflicts
    assert missing.status == "unknown" and missing.score == 0.0


def test_unresolved_condition_identity_is_visible_in_coverage() -> None:
    assessment = assess_reaction_recipe("CCBr.N>>CCN", {
        "other_components": [{"identity_status": "unresolved", "raw_identifier": "unknown additive"}],
    })
    assert assessment.coverage.condition_identity_status == "unresolved"
    assert assessment.coverage.unresolved_components == ("other_components[0]: unknown additive",)
    assert assessment.compatible and assessment.score == 0.0
