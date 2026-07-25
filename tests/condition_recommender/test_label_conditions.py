from condition_recommender.label_conditions import (
    convert_label_conditions,
    parse_process_stages,
)


def test_flat_conditions_become_registry_recipe_and_review_columns() -> None:
    converted = convert_label_conditions(
        {
            "base": "K2CO3",
            "primary_solvent": "dioxane",
            "secondary_solvent": "water",
            "ligand": "No Ligand",
            "procedure_text": "2 h at 80 °C",
        }
    )
    columns = converted.to_columns()
    payload = converted.recipe.to_dict()

    assert columns["condition_recipe_id"].startswith("RCR1:")
    assert "condition_recipe_json" not in columns
    assert columns["temperature_c"] == "80"
    assert columns["time_h"] == "2"
    assert columns["condition_identity_uncertainty"] == "false"
    assert "K2CO3 [base]" in columns["condition_display"]
    assert "No ligand" in columns["condition_display"]
    assert payload["declared_absences"] == ("ligand",)
    assert len(payload["bases"]) == 1
    assert len(payload["solvents"]) == 2
    assert all(
        component["raw_identifier"] != "No Ligand"
        for bucket in payload
        if isinstance(payload[bucket], list)
        for component in payload[bucket]
        if isinstance(component, dict) and "raw_identifier" in component
    )


def test_comma_in_registered_name_is_not_split() -> None:
    converted = convert_label_conditions(
        {"ligand": "1,10-phenanthroline", "procedure_text": ""}
    )

    assert len(converted.recipe.components) == 1
    assert converted.recipe.components[0].substance_id == "cas:66-71-7"


def test_comma_separated_components_split_only_when_each_part_resolves() -> None:
    converted = convert_label_conditions(
        {"additive": "K2CO3, NMI", "procedure_text": ""}
    )
    unresolved = convert_label_conditions(
        {"additive": "known-looking, mystery material", "procedure_text": ""}
    )

    assert len(converted.recipe.components) == 2
    assert len(unresolved.recipe.components) == 1
    assert unresolved.recipe.components[0].raw_identifier == (
        "known-looking, mystery material"
    )


def test_multi_stage_procedure_is_preserved_without_lossy_top_level_scalar() -> None:
    stages = parse_process_stages("3 h at 60 °C, then 18 h at 80 °C")
    converted = convert_label_conditions(
        {"procedure_text": "3 h at 60 °C, then 18 h at 80 °C"}
    )

    assert tuple((stage.time_h, stage.temperature_c) for stage in stages) == (
        (3.0, 60.0),
        (18.0, 80.0),
    )
    assert len(converted.recipe.stages) == 2
    assert converted.recipe.temperature_c is None
    assert converted.recipe.time_h is None


def test_recipe_identity_uses_registry_identity_not_alias_or_solvent_slot() -> None:
    first = convert_label_conditions(
        {
            "primary_solvent": "dioxane",
            "secondary_solvent": "water",
            "procedure_text": "",
        }
    )
    equivalent = convert_label_conditions(
        {
            "primary_solvent": "Water",
            "secondary_solvent": "1,4-Dioxane",
            "procedure_text": "",
        }
    )

    assert first.recipe.recipe_id == equivalent.recipe.recipe_id
