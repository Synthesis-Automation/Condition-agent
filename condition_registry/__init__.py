"""Standalone condition-substance identity and role registry."""

from .api import (
    get_registry,
    resolve_condition_components,
    resolve_identifier,
    resolve_substance,
    resolve_substance_id,
)
from .contextual_roles import resolve_contextual_component
from .curation import (
    CompoundAdditionError,
    CompoundAdditionRequest,
    CompoundAdditionResult,
    CompoundAliasInput,
    add_compound,
    update_compound,
)
from .models import (
    CONDITION_RECIPE_COMPONENT_BUCKETS,
    CONDITION_IDENTIFIER_TYPES,
    CONDITION_NAME_IDENTIFIER_TYPES,
    ConditionComponentInput,
    ConditionProcessStage,
    ContextualRoleAssignment,
    ResolutionResult,
    ResolvedConditionComponent,
    ResolvedConditionRecipe,
    RoleAssignment,
    Substance,
    SubstanceIdentifier,
)
from .recipes import (
    build_resolved_recipe,
    build_resolved_recipe_from_components,
    build_resolved_recipe_from_inputs,
)
from .resolver import ConditionRegistry
from .vocabulary import (
    ConditionDefinitionVocabulary,
    condition_registry_definition_versions,
    load_condition_vocabulary,
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
    "CONDITION_IDENTIFIER_TYPES",
    "CONDITION_NAME_IDENTIFIER_TYPES",
    "CONDITION_RECIPE_COMPONENT_BUCKETS",
    "ConditionDefinitionVocabulary",
    "ConditionRegistry",
    "CompoundAdditionError",
    "CompoundAdditionRequest",
    "CompoundAdditionResult",
    "CompoundAliasInput",
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
    "SubstanceIdentifier",
    "build_resolved_recipe",
    "build_resolved_recipe_from_components",
    "build_resolved_recipe_from_inputs",
    "condition_registry_definition_versions",
    "add_compound",
    "update_compound",
    "get_registry",
    "get_recipe_template",
    "load_recipe_template_set",
    "materialize_recipe_variant",
    "load_condition_vocabulary",
    "resolve_condition_components",
    "resolve_identifier",
    "resolve_contextual_component",
    "resolve_substance",
    "resolve_substance_id",
    "validate_registry",
    "validate_recipe_template_payload",
]
