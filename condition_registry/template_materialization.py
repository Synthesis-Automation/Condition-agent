"""Materialize explicit expert-template variants into canonical recipes."""

from __future__ import annotations

from typing import Optional

from .api import get_registry
from .models import (
    ContextualRoleAssignment,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
)
from .recipes import build_resolved_recipe_from_components
from .template_models import (
    ConditionRecipeTemplate,
    RecipeTemplateVariant,
)


def _variant(
    template: ConditionRecipeTemplate, variant_id: str
) -> RecipeTemplateVariant:
    value = next(
        (
            candidate
            for candidate in template.variants
            if candidate.variant_id == variant_id
        ),
        None,
    )
    if value is None:
        raise ValueError(
            f"Unknown recipe-template variant {template.template_id}:{variant_id}"
        )
    return value


def materialize_recipe_variant(
    template: ConditionRecipeTemplate,
    variant_id: str,
    *,
    transformation_class: Optional[str] = None,
    named_family: Optional[str] = None,
    include_draft: bool = False,
) -> ResolvedConditionRecipe:
    """Build one RCR1 recipe from an explicit, validated slot selection."""
    variant = _variant(template, variant_id)
    if not include_draft and (
        template.status != "active" or variant.status != "active"
    ):
        raise ValueError("Draft recipe-template variants require explicit opt-in")
    if template.pressure_bar is not None:
        raise ValueError(
            "ResolvedConditionRecipe does not yet represent pressure_bar"
        )
    slots = {slot.slot_id: slot for slot in template.slots}
    registry = get_registry()
    components = []
    for selection in variant.selections:
        slot = slots[selection.slot_id]
        result = registry.resolve_id(selection.substance_id)
        if result.status != "resolved" or result.substance is None:
            raise ValueError(
                f"Unresolved template substance {selection.substance_id}"
            )
        substance = result.substance
        matching_roles = tuple(
            assignment
            for assignment in substance.roles
            if assignment.role_id == slot.role_id
        )
        if not matching_roles:
            raise ValueError(
                f"Template role mismatch {selection.substance_id}:{slot.role_id}"
            )
        roles = [
            ContextualRoleAssignment(
                role_id=slot.role_id,
                confidence=1.0,
                evidence=(
                    "recipe_template_declared_role",
                    "curated_registry_identity",
                ),
            )
        ]
        roles.extend(
            ContextualRoleAssignment(
                role_id=assignment.role_id,
                confidence=0.7,
                evidence=("curated_registry_secondary_role",),
            )
            for assignment in substance.roles
            if assignment.role_id != slot.role_id
        )
        components.append(
            ResolvedConditionComponent(
                raw_identifier=selection.substance_id,
                source_field=(
                    f"recipe_template:{template.template_id}:{selection.slot_id}"
                ),
                identity_status="resolved",
                substance_id=substance.substance_id,
                canonical_name=substance.canonical_name,
                roles=tuple(roles),
                role_status="assigned",
                primary_role=slot.role_id,
                primary_role_confidence=1.0,
                cas=substance.cas,
                amount=selection.amount,
                amount_unit=selection.amount_unit,
                provenance={
                    "identity_match_kind": "exact_substance_id",
                    "recipe_template_id": template.template_id,
                    "recipe_variant_id": variant.variant_id,
                    "recipe_slot_id": selection.slot_id,
                    "transformation_class": transformation_class,
                    "named_family": named_family,
                },
            )
        )
    return build_resolved_recipe_from_components(
        components,
        temperature_c=template.temperature_c,
        time_h=template.time_h,
        concentration_m=template.concentration_m,
        atmosphere=template.atmosphere,
    )


__all__ = ["materialize_recipe_variant"]
