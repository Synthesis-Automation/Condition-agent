"""Combined source evidence and automation handoffs must remain conservative."""

import copy
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator
from referencing import Registry, Resource

from condition_registry import (
    ConditionComponentInput,
    build_resolved_recipe_from_inputs,
)
from condition_recommender.combined import (
    RecommendationSourceResult,
    build_automation_handoff,
    combine_recommendation_results,
    load_combined_recommendation_policy,
)


REACTION = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


def recipe(temperature=80.0):
    return build_resolved_recipe_from_inputs(
        (
            ConditionComponentInput(
                "K2CO3",
                source_field="base",
                identifier_type="name",
                source_role_hint="base",
                amount=2.0,
                amount_unit="equiv",
            ),
        ),
        temperature_c=temperature,
        time_h=12.0,
        atmosphere="nitrogen",
    ).to_dict()


def source(
    name="generic", *, mode="verified_signature", value=None, valid=True, score=0.1
):
    resolved = value or recipe()
    return RecommendationSourceResult(
        name,
        "ok" if valid else "abstained",
        {
            "valid": valid,
            "recommendation_mode": mode,
            "query_reaction_smiles": REACTION,
            "schema_version": "test",
            "warnings": ["SOURCE_WARNING"],
            "recommendations": [
                {
                    "recipe_id": resolved["recipe_id"],
                    "resolved_recipe": resolved,
                    "score": score,
                    "support": 2,
                    "cautions": ["Check atmosphere"],
                    "precedent_reaction_ids": ["p1"],
                }
            ],
        },
    )


def test_identical_recipes_retain_separate_source_evidence_without_adding_scores():
    combined = combine_recommendation_results(
        REACTION,
        (
            source(),
            source("weak_label", mode="weak_label_fallback", score=0.99),
        ),
    )
    assert len(combined.recommendations) == 1
    option = combined.recommendations[0]
    assert [item["source"] for item in option.evidence] == ["generic", "weak_label"]
    assert [item["recommendation"]["score"] for item in option.evidence] == [0.1, 0.99]
    assert option.evidence_kind == "verified_signature"
    assert option.cautions == ("Check atmosphere",)
    assert combined.warnings == ("SOURCE_WARNING",)


def test_different_operating_variants_are_not_merged_or_ranked_by_cross_source_score():
    combined = combine_recommendation_results(
        REACTION,
        (
            source("weak_label", value=recipe(20), score=0.99),
            source(value=recipe(80)),
        ),
    )
    assert len(combined.recommendations) == 2
    assert combined.recommendations[0].resolved_recipe["temperature_c"] == 80
    assert combined.recommendations[1].evidence_kind == "weak_label"


@pytest.mark.parametrize(
    "mode",
    [
        "reaction_core_review",
        "ambiguous_edit_hypotheses",
        "structure_fallback",
        "unknown",
    ],
)
def test_ambiguous_and_review_modes_never_receive_verified_evidence_badge(mode):
    combined = combine_recommendation_results(REACTION, (source(mode=mode),))
    assert combined.recommendations[0].evidence_kind == "structure_review"


def test_conflicting_abstained_output_cannot_leak_stale_recommendations():
    combined = combine_recommendation_results(REACTION, (source(valid=False),))
    assert not combined.valid
    assert not combined.recommendations
    assert combined.sources[0].result["recommendations"]  # kept only as audit


def test_partial_source_failure_retains_other_evidence_and_status():
    combined = combine_recommendation_results(
        REACTION,
        (
            RecommendationSourceResult(
                "generic", "unavailable", {}, "Index unavailable"
            ),
            source("weak_label", mode="weak_label_fallback"),
        ),
    )
    assert combined.valid
    assert combined.recommendations[0].evidence_kind == "weak_label"
    assert combined.sources[0].status == "unavailable"


def test_quantities_stages_and_definition_versions_survive_automation_export():
    resolved = recipe()
    resolved["stages"] = [
        {"stage_index": 1, "temperature_c": 20, "time_h": 1},
        {"stage_index": 2, "temperature_c": 80, "time_h": 12},
    ]
    combined = combine_recommendation_results(REACTION, (source(value=resolved),))
    before = copy.deepcopy(combined.to_dict())
    ids = (combined.recommendations[0].option_id,)
    handoff = build_automation_handoff(combined, ids)
    assert handoff == build_automation_handoff(combined, ids)
    assert combined.to_dict() == before
    assert handoff["execution_ready"] is False
    assert handoff["robot_target"] is None
    experiment = handoff["experiments"][0]
    assert experiment["resolved_recipe"] == resolved
    protocol = experiment["protocol"]
    assert [item["temperature_c"] for item in protocol["operations"]] == [20, 80]
    condition = next(
        item for item in protocol["materials"] if item["category"] == "condition"
    )
    assert (condition["amount"], condition["amount_unit"]) == (2, "equiv")
    assert "ordered_operations" in protocol["missing_required_fields"]
    assert "materials.reactant_001.amount" in protocol["missing_required_fields"]
    root = Path(__file__).resolve().parents[2]
    schema = json.loads(
        (
            root / "condition_recommender/definitions/automation_handoff.v1.schema.json"
        ).read_text()
    )
    protocol_schema = json.loads(
        (
            root / "condition_registry/definitions/synthesis_protocol.v1.schema.json"
        ).read_text()
    )
    Draft202012Validator.check_schema(schema)
    registry = Registry().with_resource(
        protocol_schema["$id"], Resource.from_contents(protocol_schema)
    )
    Draft202012Validator(schema, registry=registry).validate(
        json.loads(json.dumps(handoff))
    )


def test_export_uses_confirmed_effective_query_and_keeps_original_query():
    original = source()
    result = dict(original.result)
    result["effective_query_reaction_smiles"] = REACTION + ".O"
    combined = combine_recommendation_results(
        REACTION, (RecommendationSourceResult("generic", "ok", result),)
    )
    assert combined.query_reaction_smiles == REACTION
    assert (
        combined.recommendations[0].synthesis_protocol["reaction_smiles"]
        == REACTION + ".O"
    )


def test_export_rejects_empty_unknown_and_duplicate_selection():
    combined = combine_recommendation_results(REACTION, (source(),))
    option_id = combined.recommendations[0].option_id
    for selection in ((), ("unknown",), (option_id, option_id)):
        with pytest.raises(ValueError):
            build_automation_handoff(combined, selection)


def test_policy_and_option_ids_are_deterministic():
    assert load_combined_recommendation_policy()["default_shortlist_size"] == 3
    assert (
        combine_recommendation_results(REACTION, (source(),)).to_dict()
        == combine_recommendation_results(REACTION, (source(),)).to_dict()
    )
