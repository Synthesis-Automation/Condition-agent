"""Standalone condition-substance identity and role registry."""

from .api import get_registry, resolve_condition_components, resolve_substance
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
from .validation import validate_registry

__all__ = [
    "ContextualRoleAssignment",
    "ResolutionResult",
    "ResolvedConditionComponent",
    "ResolvedConditionRecipe",
    "RoleAssignment",
    "Substance",
    "build_resolved_recipe",
    "get_registry",
    "resolve_condition_components",
    "resolve_contextual_component",
    "resolve_substance",
    "validate_registry",
]
