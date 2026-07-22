"""Standalone condition-substance identity and role registry."""

from .api import (
    get_registry,
    resolve_condition_components,
    resolve_substance,
    resolve_substance_id,
)
from .contextual_roles import resolve_contextual_component
from .models import (
    ContextualRoleAssignment,
    ResolutionResult,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
    RoleAssignment,
    Substance,
)
from .recipes import build_resolved_recipe
from .template_loader import (
    get_recipe_template,
    load_recipe_template_set,
    validate_recipe_template_payload,
)
from .template_models import (
    ConditionRecipeTemplate,
    ConditionRecipeTemplateSet,
    RecipeTemplateOption,
    RecipeTemplateSlot,
)
from .validation import validate_registry

__all__ = [
    "ContextualRoleAssignment",
    "ConditionRecipeTemplate",
    "ConditionRecipeTemplateSet",
    "ResolutionResult",
    "ResolvedConditionComponent",
    "ResolvedConditionRecipe",
    "RoleAssignment",
    "RecipeTemplateOption",
    "RecipeTemplateSlot",
    "Substance",
    "build_resolved_recipe",
    "get_registry",
    "get_recipe_template",
    "load_recipe_template_set",
    "resolve_condition_components",
    "resolve_contextual_component",
    "resolve_substance",
    "resolve_substance_id",
    "validate_registry",
    "validate_recipe_template_payload",
]
