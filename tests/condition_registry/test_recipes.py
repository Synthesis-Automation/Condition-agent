from condition_registry import (
    ConditionComponentInput,
    ConditionProcessStage,
    build_resolved_recipe,
    build_resolved_recipe_from_inputs,
    resolve_contextual_component,
)
from condition_registry.contextual_roles import load_role_resolution_rules
from condition_registry.loader import load_role_definitions


def test_multi_role_substance_uses_source_and_reaction_role_context() -> None:
    reagent = resolve_contextual_component(
        "121-44-8",
        source_field="reagent_cas",
        preferred_roles=("metal_catalyst", "ligand", "base"),
    )
    solvent = resolve_contextual_component(
        "121-44-8",
        source_field="solvent_cas",
        preferred_roles=("metal_catalyst", "ligand", "base"),
    )

    assert reagent.primary_role == "base"
    assert solvent.primary_role == "solvent"
    assert {role.role_id for role in reagent.roles} == {"base", "solvent"}
    assert reagent.role_status == "assigned"


def test_resolved_recipe_groups_components_by_contextual_role() -> None:
    recipe = build_resolved_recipe(
        {
            "catalyst_cas": ("14221-01-3",),
            "reagent_cas": ("584-08-7",),
            "solvent_cas": ("7732-18-5",),
        },
        preferred_roles=("metal_catalyst", "ligand", "base"),
        temperature_c=80.0,
        time_h=12.0,
    )

    assert recipe.recipe_id.startswith("RCR2:")
    assert recipe.recipe_core_id.startswith("RCORE2:")
    assert recipe.catalysts[0].primary_role == "metal_catalyst"
    assert recipe.bases[0].primary_role == "base"
    assert recipe.solvents[0].primary_role == "solvent"
    assert recipe.solvents[0].cas == "7732-18-5"
    assert len(recipe.components) == 3
    assert recipe.warnings == ()


def test_recipe_identity_is_input_order_invariant() -> None:
    first = build_resolved_recipe(
        {
            "reagent_cas": ("584-08-7", "121-44-8"),
            "solvent_cas": ("7732-18-5",),
        },
        preferred_roles=("metal_catalyst", "ligand", "base"),
    )
    reordered = build_resolved_recipe(
        {
            "solvent_cas": ("7732-18-5",),
            "reagent_cas": ("121-44-8", "584-08-7"),
        },
        preferred_roles=("metal_catalyst", "ligand", "base"),
    )

    assert first.recipe_id == reordered.recipe_id


def test_unresolved_identity_is_retained_in_other_components() -> None:
    recipe = build_resolved_recipe(
        {"catalyst_cas": ("7440-06-4",)},
        preferred_roles=("metal_catalyst", "ligand", "base"),
    )

    assert len(recipe.other_components) == 1
    component = recipe.other_components[0]
    assert component.raw_identifier == "7440-06-4"
    assert component.identity_status == "unresolved"
    assert component.role_status == "unassigned"
    assert component.primary_role is None
    assert component.cas is None
    assert "CONDITION_IDENTITY_UNCERTAINTY" in recipe.warnings


def test_operating_conditions_are_part_of_recipe_identity() -> None:
    ambient = build_resolved_recipe(
        {"solvent_cas": ("7732-18-5",)}, temperature_c=25.0
    )
    heated = build_resolved_recipe(
        {"solvent_cas": ("7732-18-5",)}, temperature_c=80.0
    )

    assert ambient.recipe_id != heated.recipe_id
    assert ambient.recipe_core_id == heated.recipe_core_id


def test_component_amount_changes_variant_but_not_recipe_core() -> None:
    first = build_resolved_recipe_from_inputs(
        (
            ConditionComponentInput(
                "K2CO3",
                source_field="base",
                identifier_type="name",
                amount=1.0,
                amount_unit="equiv",
            ),
        )
    )
    second = build_resolved_recipe_from_inputs(
        (
            ConditionComponentInput(
                "K2CO3",
                source_field="base",
                identifier_type="name",
                amount=2.0,
                amount_unit="equiv",
            ),
        )
    )

    assert first.recipe_id != second.recipe_id
    assert first.recipe_core_id == second.recipe_core_id


def test_duplicate_identity_across_source_fields_is_not_double_counted() -> None:
    recipe = build_resolved_recipe(
        {
            "reagent_cas": ("7732-18-5",),
            "solvent_cas": ("7732-18-5",),
        }
    )

    assert len(recipe.components) == 1
    assert recipe.components[0].primary_role == "solvent"
    assert recipe.components[0].provenance["source_fields"] == (
        "reagent_cas",
        "solvent_cas",
    )
    assert "DUPLICATE_SOURCE_IDENTITY_MERGED" in recipe.warnings


def test_role_resolution_definition_references_known_roles_and_buckets() -> None:
    taxonomy_roles = {item["id"] for item in load_role_definitions()["roles"]}
    rules = load_role_resolution_rules()

    assert set(rules["role_buckets"]) <= taxonomy_roles
    assert set(rules["source_role_fallbacks"].values()) <= taxonomy_roles
    assert set(rules["role_buckets"].values()) == {
        "catalysts",
        "ligands",
        "bases",
        "acids",
        "condensation_agents",
        "oxidants",
        "reductants",
        "additives",
        "solvents",
    }


def test_contextual_component_resolves_name_and_retains_role_hint() -> None:
    component = resolve_contextual_component(
        "K2CO3",
        source_field="base",
        identifier_type="name",
        source_role_hint="base",
    )

    assert component.substance_id == "cas:584-08-7"
    assert component.canonical_name == "Potassium carbonate"
    assert component.primary_role == "base"
    assert component.source_role_hint == "base"
    assert component.provenance["identifier_type"] == "name"
    assert component.provenance["identity_identifier_type"] == "abbreviation"
    assert "identity_identifier_source" not in component.provenance


def test_source_role_hint_does_not_override_registry_chemistry() -> None:
    component = resolve_contextual_component(
        "water",
        source_field="base",
        identifier_type="name",
        source_role_hint="base",
    )

    assert component.role_status == "conflicting"
    assert component.primary_role is None
    assert "SOURCE_ROLE_CONFLICT" in component.warnings


def test_unresolved_name_uses_explicit_hint_fallback_with_uncertainty() -> None:
    component = resolve_contextual_component(
        "unknown ligand",
        source_field="ligand",
        identifier_type="name",
        source_role_hint="ligand",
    )

    assert component.identity_status == "unresolved"
    assert component.primary_role == "ligand"
    assert component.primary_role_confidence < 0.7
    assert component.roles[0].evidence[0] == "source_role_hint"


def test_recipe_identity_includes_stages_and_declared_absences() -> None:
    inputs = (
        ConditionComponentInput(
            "water",
            source_field="primary_solvent",
            identifier_type="name",
            source_role_hint="solvent",
        ),
    )
    one_stage = build_resolved_recipe_from_inputs(
        inputs,
        stages=(ConditionProcessStage(1, temperature_c=25.0, time_h=1.0),),
        declared_absences=("ligand",),
    )
    two_stage = build_resolved_recipe_from_inputs(
        inputs,
        stages=(
            ConditionProcessStage(1, temperature_c=25.0, time_h=1.0),
            ConditionProcessStage(2, temperature_c=80.0, time_h=2.0),
        ),
        declared_absences=("ligand",),
    )
    ligand_present_or_unspecified = build_resolved_recipe_from_inputs(
        inputs,
        stages=(ConditionProcessStage(1, temperature_c=25.0, time_h=1.0),),
    )

    assert one_stage.recipe_id != two_stage.recipe_id
    assert one_stage.recipe_id != ligand_present_or_unspecified.recipe_id
    assert one_stage.declared_absences == ("ligand",)
