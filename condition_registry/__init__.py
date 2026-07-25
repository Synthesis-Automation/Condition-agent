"""Standalone condition-substance identity and role registry."""

from .api import (
    get_registry,
    resolve_condition_components,
    resolve_substance,
    resolve_substance_id,
)
from .contextual_roles import resolve_contextual_component
from .models import (
    ConditionComponentInput,
    ConditionProcessStage,
    ContextualRoleAssignment,
    ResolutionResult,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
    RoleAssignment,
    Substance,
)
from .recipes import (
    build_resolved_recipe,
    build_resolved_recipe_from_components,
    build_resolved_recipe_from_inputs,
)
from .template_loader import (
    get_recipe_template,
    load_recipe_template_set,
    validate_recipe_template_payload,
)
from .template_materialization import materialize_recipe_variant
from .template_models import (
    ConditionRecipeTemplate,
    ConditionRecipeTemplateSet,
    RecipeTemplateOption,
    RecipeTemplatePartnerAmount,
    RecipeTemplateSelection,
    RecipeTemplateSlot,
    RecipeTemplateVariant,
)
from .validation import validate_registry

__all__ = [
    "ContextualRoleAssignment",
    "ConditionComponentInput",
    "ConditionProcessStage",
    "ConditionRecipeTemplate",
    "ConditionRecipeTemplateSet",
    "ResolutionResult",
    "ResolvedConditionComponent",
    "ResolvedConditionRecipe",
    "RoleAssignment",
    "RecipeTemplateOption",
    "RecipeTemplatePartnerAmount",
    "RecipeTemplateSelection",
    "RecipeTemplateSlot",
    "RecipeTemplateVariant",
    "Substance",
    "build_resolved_recipe",
    "build_resolved_recipe_from_components",
    "build_resolved_recipe_from_inputs",
    "get_registry",
    "get_recipe_template",
    "load_recipe_template_set",
    "materialize_recipe_variant",
    "resolve_condition_components",
    "resolve_contextual_component",
    "resolve_substance",
    "resolve_substance_id",
    "validate_registry",
    "validate_recipe_template_payload",
]
