"""General condition-source capability matching regressions."""

from dataclasses import asdict

from condition_registry import (
    ConditionComponentInput,
    build_resolved_recipe,
    build_resolved_recipe_from_inputs,
)
from reactive_taxonomy import featurize_reaction

from condition_recommender.fragment_source_support import (
    assess_fragment_source_support,
    fragment_source_support_is_complete,
    load_fragment_source_capabilities,
)


def _requirements(reaction_smiles: str):
    descriptor = featurize_reaction(reaction_smiles).fallback_descriptor
    assert descriptor is not None
    return descriptor.source_requirements


def test_fragment_source_capability_definition_is_valid() -> None:
    rules = load_fragment_source_capabilities()

    assert rules["definition_id"] == "fragment_source_capabilities.v1"
    assert len(rules["capabilities"]) >= 6


def test_single_atom_source_requires_curated_transfer_capability() -> None:
    requirements = _requirements("CC(C)O>>CC(C)F")
    capable = build_resolved_recipe(
        {"reagent_cas": ("429-41-4",)},
    )
    merely_fluorinated = build_resolved_recipe(
        {"reagent_cas": ("920-66-1",)},
    )

    supported = assess_fragment_source_support(requirements, capable)
    unsupported = assess_fragment_source_support(
        requirements,
        merely_fluorinated,
    )

    assert supported[0].status == "supported"
    assert supported[0].capability_ids == ("single_atom_fluorine_source.v1",)
    assert fragment_source_support_is_complete(requirements, supported)
    assert unsupported[0].status == "unsupported"
    assert not fragment_source_support_is_complete(requirements, unsupported)


def test_multi_atom_cyanide_requirement_uses_same_general_contract() -> None:
    requirements = _requirements("Brc1ccccc1>>N#Cc1ccccc1")
    recipe = build_resolved_recipe(
        {"reagent_cas": ("143-33-9",)},
    )

    support = assess_fragment_source_support(requirements, recipe)

    requirement = asdict(requirements[0])
    assert requirement["rooted_fragment_smiles"] == "*C#N"
    assert requirement["element_counts"] == {"C": 1, "N": 1}
    assert support[0].status == "supported"
    assert support[0].capability_ids == ("cyanide_fragment_source.v1",)


def test_structured_reagent_field_supports_azide_requirement() -> None:
    requirements = _requirements(
        "Fc1ccccc1>>[N-]=[N+]=Nc1ccccc1"
    )
    recipe = build_resolved_recipe_from_inputs(
        (
            ConditionComponentInput(
                raw_identifier="26628-22-8",
                source_field="reagents_json",
                identifier_type="cas",
                source_role_hint="reagent",
            ),
        ),
    )

    support = assess_fragment_source_support(requirements, recipe)

    assert support[0].status == "supported"
    assert support[0].component_substance_ids == ("cas:26628-22-8",)
    assert support[0].capability_ids == ("azide_fragment_source.v1",)
