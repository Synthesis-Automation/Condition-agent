from types import SimpleNamespace

from condition_recommender.compatibility import (
    assess_recipe_compatibility,
    filter_compatible_precedents,
    load_compatibility_rules,
)


def _signature(*tags: str) -> dict:
    return {
        "spectator_groups": [
            {"group_id": "test_group", "tags": list(tags)}
        ]
    }


def _component(family_id: str, role_id: str) -> dict:
    return {
        "roles": [{"family_id": family_id, "role_id": role_id}]
    }


def _recipe(**buckets: list[dict]) -> dict:
    return buckets


def test_oxidation_sensitive_spectator_rejects_oxidant_recipe() -> None:
    assessment = assess_recipe_compatibility(
        _signature("oxidation_sensitive"),
        _recipe(oxidants=[_component("peroxides_hydroperoxides", "oxidant")]),
    )
    assert not assessment.compatible
    assert assessment.score == 0.0
    assert assessment.hard_conflicts == ("oxidation_sensitive_with_oxidant",)


def test_moisture_sensitive_spectator_rejects_water_recipe() -> None:
    assessment = assess_recipe_compatibility(
        _signature("moisture_sensitive"),
        _recipe(solvents=[_component("water", "solvent")]),
    )
    assert not assessment.compatible
    assert assessment.hard_conflicts == ("moisture_sensitive_with_water",)


def test_overlapping_coordination_risks_use_only_strongest_penalty() -> None:
    assessment = assess_recipe_compatibility(
        _signature(
            "metal_binding", "strong_metal_binding", "palladium_poisoning"
        ),
        _recipe(catalysts=[_component("pd_zero_sources", "metal_catalyst")]),
    )
    assert assessment.compatible
    assert assessment.score == 0.75
    assert assessment.penalty_ids == ("palladium_poisoning_with_pd_catalyst",)
    assert assessment.evidence == (
        "A spectator group may poison the palladium catalyst",
    )


def test_hot_water_uses_stronger_hydrolysis_penalty() -> None:
    query = _signature("hydrolysis_sensitive")
    water = [_component("water", "solvent")]
    ambient = assess_recipe_compatibility(
        query, {"solvents": water, "temperature_c": 25.0}
    )
    hot = assess_recipe_compatibility(
        query, {"solvents": water, "temperature_c": 80.0}
    )
    assert ambient.score == 0.75
    assert ambient.penalty_ids == ("hydrolysis_sensitive_with_water",)
    assert hot.score == 0.6
    assert hot.penalty_ids == ("hydrolysis_sensitive_with_hot_water",)


def test_only_unchanged_spectator_tags_drive_compatibility() -> None:
    signature = {
        "spectator_groups": [],
        "partners": [
            {"nearby_groups": [{"tags": ["oxidation_sensitive"]}]}
        ],
    }
    assessment = assess_recipe_compatibility(
        signature,
        _recipe(oxidants=[_component("peroxides_hydroperoxides", "oxidant")]),
    )
    assert assessment.compatible
    assert assessment.score == 1.0


def test_filter_retains_exclusion_evidence() -> None:
    safe = SimpleNamespace(resolved_recipe=_recipe(bases=[]))
    unsafe = SimpleNamespace(
        resolved_recipe=_recipe(
            reductants=[_component("borohydrides", "reductant")]
        )
    )
    accepted, excluded = filter_compatible_precedents(
        _signature("reduction_sensitive"), (safe, unsafe)
    )
    assert tuple(item[0] for item in accepted) == (safe,)
    assert tuple(item[0] for item in excluded) == (unsafe,)
    assert excluded[0][1].hard_conflicts == (
        "reduction_sensitive_with_reductant",
    )


def test_oxidative_atmosphere_conflict_does_not_match_anaerobic() -> None:
    query = _signature("oxidation_sensitive")
    aerobic = assess_recipe_compatibility(query, {"atmosphere": "aerobic"})
    anaerobic = assess_recipe_compatibility(query, {"atmosphere": "anaerobic"})
    assert not aerobic.compatible
    assert aerobic.hard_conflicts == (
        "oxidation_sensitive_in_oxidative_atmosphere",
    )
    assert anaerobic.compatible


def test_mandatory_catalyst_is_hard_unless_identity_is_unresolved() -> None:
    query = {
        **_signature(),
        "named_family": "chan_lam",
        "family_confidence": 1.0,
    }
    missing = assess_recipe_compatibility(query, {})
    assert not missing.compatible
    assert missing.hard_conflicts == ("metal_coupling_requires_catalyst",)

    unresolved = assess_recipe_compatibility(
        query,
        {
            "other_components": [
                {
                    "identity_status": "unresolved",
                    "roles": [],
                }
            ]
        },
    )
    assert unresolved.compatible
    assert unresolved.score == 0.88
    assert unresolved.penalty_ids == ("metal_coupling_requires_catalyst",)

    low_confidence = assess_recipe_compatibility(
        {**query, "family_confidence": 0.5}, {}
    )
    assert low_confidence.compatible
    assert low_confidence.score == 1.0


def test_compatibility_definition_has_unique_rule_ids() -> None:
    rules = load_compatibility_rules()
    all_rules = (
        rules["hard_conflicts"]
        + rules["soft_penalties"]
        + rules["regime_requirements"]
    )
    ids = [rule["id"] for rule in all_rules]
    assert len(ids) == len(set(ids))
    assert all(0.0 < rule["penalty"] < 1.0 for rule in rules["soft_penalties"])
