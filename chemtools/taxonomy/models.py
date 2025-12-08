"""
Dataclasses describing the unified taxonomy entities.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass(slots=True)
class TaxonomyManifest:
    """Metadata stored alongside the taxonomy payload."""

    schema_version: str
    taxonomy_version: str
    generated_at: str
    source: Optional[str] = None
    notes: Optional[str] = None


@dataclass(slots=True)
class ReactionCategory:
    id: str
    name: str
    description: Optional[str] = None


@dataclass(slots=True)
class ReactionReactantRequirement:
    reactant_type_id: str
    stoichiometry: Optional[str] = None
    notes: Optional[str] = None
    original_tokens: List[str] = field(default_factory=list)


@dataclass(slots=True)
class ReactionRoleRequirement:
    role_id: str
    required: bool = True
    default_family_id: Optional[str] = None
    notes: Optional[str] = None


@dataclass(slots=True)
class ReactionPattern:
    """A SMARTS pattern for reaction matching."""
    id: str
    smarts: str
    source: Optional[str] = None
    scope: Optional[str] = None


@dataclass(slots=True)
class ReactionExample:
    """Example reactants for a reaction type."""
    reactant1: str
    reactant2: Optional[str] = None
    source: Optional[str] = None


@dataclass(slots=True)
class ReactionType:
    id: str
    category_id: str
    name: str
    description: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    reactants: List[ReactionReactantRequirement] = field(default_factory=list)
    required_roles: List[ReactionRoleRequirement] = field(default_factory=list)
    patterns: List[ReactionPattern] = field(default_factory=list)
    examples: List[ReactionExample] = field(default_factory=list)
    catalysts: List[str] = field(default_factory=list)
    conditions: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    source_ids: List[str] = field(default_factory=list)


@dataclass(slots=True)
class ReactantTypeMember:
    id: str
    name: str
    smarts: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ReactantType:
    id: str
    name: str
    description: Optional[str] = None
    category: Optional[str] = None
    smarts: Optional[str] = None
    members: List[ReactantTypeMember] = field(default_factory=list)
    aliases: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ReagentRole:
    id: str
    name: str
    priority: int = 100
    default_family_id: Optional[str] = None
    description: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ReagentFamily:
    id: str
    role_id: str
    name: str
    description: Optional[str] = None
    precedence: Optional[int] = None
    include_smarts: List[str] = field(default_factory=list)
    exclude_smarts: List[str] = field(default_factory=list)
    required_props: Dict[str, Any] = field(default_factory=dict)
    keywords: List[str] = field(default_factory=list)
    examples_pos: List[str] = field(default_factory=list)
    examples_neg: List[str] = field(default_factory=list)
    notes: Optional[str] = None


@dataclass(slots=True)
class AliasRecord:
    alias: str
    entity_type: str
    entity_id: str
    source: Optional[str] = None
    notes: Optional[str] = None
