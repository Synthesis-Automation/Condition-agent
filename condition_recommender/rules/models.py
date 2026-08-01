"""Typed contracts for declarative expert-condition recommendation rules."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Optional, Tuple


RuleStatus = Literal["draft", "active", "retired"]
RuleTier = Literal["specific", "fallback"]
RuleKind = Literal["default", "structural_override", "fallback"]
RecipeVariantStatus = Literal["draft", "active", "retired"]


@dataclass(frozen=True)
class PartnerConstraint:
    """Allowlisted predicates for one interpreted reaction partner."""

    role: str
    site_type: str
    anchor_contexts_any: Tuple[str, ...] = ()
    handle_tokens_any: Tuple[str, ...] = ()
    center_tokens_any: Tuple[str, ...] = ()
    derived_families_any: Tuple[str, ...] = ()
    h_count_min: Optional[int] = None
    h_count_max: Optional[int] = None
    retained_contexts_any: Tuple[str, ...] = ()
    retained_contexts_all: Tuple[str, ...] = ()
    retained_contexts_allowed: Tuple[str, ...] = ()
    availability_any: Tuple[str, ...] = ()
    ortho_occupancy_min: Optional[int] = None
    ortho_occupancy_max: Optional[int] = None
    alpha_branched_group_count_min: Optional[int] = None
    alpha_branched_group_count_max: Optional[int] = None
    context_kinds_any: Tuple[str, ...] = ()
    ring_families_any: Tuple[str, ...] = ()
    steric_accessibility_any: Tuple[str, ...] = ()
    ortho_burden_classes_any: Tuple[str, ...] = ()
    electronic_axes_any: Tuple[str, ...] = ()
    alkyl_substitutions_any: Tuple[str, ...] = ()
    beta_hydrogen_statuses_any: Tuple[str, ...] = ()
    lone_pair_availability_any: Tuple[str, ...] = ()
    reactivity_modifiers_any: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RuleScope:
    """Reaction-event scope required before partner predicates are evaluated."""

    transformation_classes_any: Tuple[str, ...]
    event_scopes_any: Tuple[str, ...] = ("single_event",)
    evidence_qualities_any: Tuple[str, ...] = (
        "exact_product_reconstruction",
    )
    reaction_scopes_any: Tuple[str, ...] = ("intermolecular",)


@dataclass(frozen=True)
class RuleSelection:
    """Deterministic specificity group and priority."""

    group: str
    tier: RuleTier
    priority: int


@dataclass(frozen=True)
class RuleRecommendation:
    """Reference from an applicability rule to a registry-owned recipe template."""

    recipe_template_id: str


@dataclass(frozen=True)
class ConditionRule:
    """One declarative structural rule with no executable chemistry payload."""

    rule_id: str
    status: RuleStatus
    rule_kind: RuleKind
    scope: RuleScope
    partner_constraints: Tuple[PartnerConstraint, ...]
    selection: RuleSelection
    recommendations: Tuple[RuleRecommendation, ...]
    rationale: str
    cautions: Tuple[str, ...] = ()
    provenance: Dict[str, object] = None  # type: ignore[assignment]
    schema_version: str = "2.0"

    def __post_init__(self) -> None:
        if self.provenance is None:
            object.__setattr__(self, "provenance", {})

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ConditionRuleSet:
    """Versioned rules plus their deterministic tier-selection policy."""

    definition_id: str
    selection_mode: Literal["first_nonempty_tier_highest_priority"]
    tier_order: Tuple[RuleTier, ...]
    rules: Tuple[ConditionRule, ...]
    schema_version: str = "2.0"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class PartnerRuleFacts:
    """Stable structural facts projected from one selected taxonomy role."""

    role: str
    component_index: int
    site_id: str
    site_type: str
    availability: str
    anchor_context: Optional[str]
    handle_token: Optional[str]
    center_token: Optional[str]
    derived_family: Optional[str]
    h_count: Optional[int]
    retained_contexts: Tuple[str, ...]
    ortho_occupancy: Optional[int] = None
    alpha_branched_group_count: int = 0
    environment_flags: Tuple[str, ...] = ()
    context_kind: Optional[str] = None
    ring_family: Optional[str] = None
    steric_accessibility: Optional[str] = None
    ortho_burden_class: Optional[str] = None
    electronic_axis: Optional[str] = None
    alkyl_substitution: Optional[str] = None
    beta_hydrogen_status: Optional[str] = None
    lone_pair_availability: Optional[str] = None
    reactivity_modifiers: Tuple[str, ...] = ()


@dataclass(frozen=True)
class RuleQueryFacts:
    """Verified event facts consumed by condition rules."""

    signature_id: str
    reaction_signature_schema_version: str
    transformation_class: str
    event_scope: str
    evidence_quality: str
    reaction_scope: str
    partners: Tuple[PartnerRuleFacts, ...]
    taxonomy_definition_versions: Tuple[Tuple[str, str], ...]


@dataclass(frozen=True)
class RuleMatchTrace:
    """Complete structural match result for one rule."""

    rule_id: str
    rule_status: RuleStatus
    rule_kind: RuleKind
    matched: bool
    eligible: bool
    selection_group: str
    selection_tier: RuleTier
    priority: int
    evidence: Tuple[str, ...] = ()
    failures: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RuleConditionRecommendation:
    """One selected expert rule and its registry-owned recipe template."""

    rank: int
    rule_id: str
    rule_status: RuleStatus
    rule_kind: RuleKind
    recipe_template_id: str
    recipe_template: Dict[str, Any]
    compatible_variants: Tuple["RuleRecipeVariantAssessment", ...]
    selection_group: str
    selection_tier: RuleTier
    rationale: str
    match_evidence: Tuple[str, ...]
    cautions: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RuleRecipeVariantAssessment:
    """Compatibility result for one explicit materialized recipe variant."""

    rule_id: str
    recipe_template_id: str
    variant_id: str
    variant_status: RecipeVariantStatus
    rank: int
    definition_priority: int
    recipe_id: str
    resolved_recipe: Dict[str, Any]
    compatible: bool
    compatibility_score: float
    hard_conflicts: Tuple[str, ...] = ()
    penalty_ids: Tuple[str, ...] = ()
    compatibility_evidence: Tuple[str, ...] = ()
    cautions: Tuple[str, ...] = ()
    compatibility_definition_id: str = "compatibility.v1"
    schema_version: str = "1.2"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RuleRecommendationResult:
    """Auditable expert-rule result for one verified reaction event."""

    query_reaction_smiles: str
    valid: bool
    query_signature_id: Optional[str] = None
    reaction_signature_schema_version: Optional[str] = None
    transformation_class: Optional[str] = None
    taxonomy_definition_versions: Tuple[Tuple[str, str], ...] = ()
    selected_tiers: Tuple[Tuple[str, str], ...] = ()
    recommendations: Tuple[RuleConditionRecommendation, ...] = ()
    excluded_variants: Tuple[RuleRecipeVariantAssessment, ...] = ()
    match_traces: Tuple[RuleMatchTrace, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    rule_definition_id: Optional[str] = None
    rule_definition_schema_version: Optional[str] = None
    recipe_template_definition_id: Optional[str] = None
    recipe_template_schema_version: Optional[str] = None
    schema_version: str = "2.0"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


__all__ = [
    "ConditionRule",
    "ConditionRuleSet",
    "PartnerConstraint",
    "PartnerRuleFacts",
    "RuleConditionRecommendation",
    "RuleKind",
    "RuleMatchTrace",
    "RuleQueryFacts",
    "RuleRecipeVariantAssessment",
    "RuleRecommendationResult",
    "RecipeVariantStatus",
    "RuleRecommendation",
    "RuleScope",
    "RuleSelection",
    "RuleStatus",
    "RuleTier",
]
