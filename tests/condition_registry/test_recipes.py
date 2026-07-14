from condition_registry import build_resolved_recipe, resolve_contextual_component
from condition_registry.contextual_roles import load_role_resolution_rules
from condition_registry.loader import load_taxonomy


def test_multi_role_substance_uses_source_and_transformation_context() -> None:
    reagent = resolve_contextual_component(
        "121-44-8",
        source_field="reagent_cas",
        transformation_class="c_c_transfer_coupling",
    )
    solvent = resolve_contextual_component(
        "121-44-8",
        source_field="solvent_cas",
        transformation_class="c_c_transfer_coupling",
    )

    assert reagent.primary_role == "base"
    assert solvent.primary_role == "solvent"
    assert {role.role_id for role in reagent.roles} == {"base", "solvent"}
    assert "MULTIPLE_POSSIBLE_ROLES" in reagent.warnings


def test_resolved_recipe_groups_components_by_contextual_role() -> None:
    recipe = build_resolved_recipe(
        {
            "catalyst_cas": ("14221-01-3",),
            "reagent_cas": ("584-08-7",),
            "solvent_cas": ("7732-18-5",),
        },
        transformation_class="c_c_transfer_coupling",
        named_family="suzuki_miyaura",
        temperature_c=80.0,
        time_h=12.0,
    )

    assert recipe.recipe_id.startswith("RCR1:")
    assert recipe.catalysts[0].primary_role == "metal_catalyst"
    assert recipe.bases[0].primary_role == "base"
    assert recipe.solvents[0].primary_role == "solvent"
    assert len(recipe.components) == 3
    assert recipe.warnings == ()


def test_recipe_identity_is_input_order_invariant() -> None:
    first = build_resolved_recipe(
        {
            "reagent_cas": ("584-08-7", "121-44-8"),
            "solvent_cas": ("7732-18-5",),
        },
        transformation_class="c_c_transfer_coupling",
    )
    reordered = build_resolved_recipe(
        {
            "solvent_cas": ("7732-18-5",),
            "reagent_cas": ("121-44-8", "584-08-7"),
        },
        transformation_class="c_c_transfer_coupling",
    )

    assert first.recipe_id == reordered.recipe_id


def test_unresolved_identity_is_retained_in_other_components() -> None:
    recipe = build_resolved_recipe(
        {"catalyst_cas": ("7440-06-4",)},
        transformation_class="c_c_transfer_coupling",
    )

    assert len(recipe.other_components) == 1
    component = recipe.other_components[0]
    assert component.raw_identifier == "7440-06-4"
    assert component.identity_status == "unresolved"
    assert component.primary_role == "other_reagent"
    assert "CONDITION_IDENTITY_UNCERTAINTY" in recipe.warnings


def test_operating_conditions_are_part_of_recipe_identity() -> None:
    ambient = build_resolved_recipe(
        {"solvent_cas": ("7732-18-5",)}, temperature_c=25.0
    )
    heated = build_resolved_recipe(
        {"solvent_cas": ("7732-18-5",)}, temperature_c=80.0
    )

    assert ambient.recipe_id != heated.recipe_id


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
    taxonomy_roles = {item["id"] for item in load_taxonomy()["roles"]}
    rules = load_role_resolution_rules()

    assert set(rules["role_buckets"]) <= taxonomy_roles
    assert set(rules["source_fallback_roles"].values()) <= taxonomy_roles
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
        "other_components",
    }
