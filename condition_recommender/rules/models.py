"""Typed contracts for declarative expert-condition recommendation rules."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Literal, Optional, Tuple


RuleStatus = Literal["draft", "active", "retired"]
RuleTier = Literal["specific", "fallback"]


@dataclass(frozen=True)
class PartnerConstraint:
    """Allowlisted structural predicates for one grammar-assigned partner."""

    role: str
    site_type: str
    anchor_contexts_any: Tuple[str, ...] = ()
    center_tokens_any: Tuple[str, ...] = ()
    derived_families_any: Tuple[str, ...] = ()
    h_count_min: Optional[int] = None
    h_count_max: Optional[int] = None
    retained_contexts_all: Tuple[str, ...] = ()
    retained_contexts_allowed: Tuple[str, ...] = ()
    availability_any: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RuleScope:
    """Reaction-event scope required before partner predicates are evaluated."""

    transformation_classes_any: Tuple[str, ...]
    event_scopes_any: Tuple[str, ...] = ("single_event",)


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
    scope: RuleScope
    partner_constraints: Tuple[PartnerConstraint, ...]
    selection: RuleSelection
    recommendations: Tuple[RuleRecommendation, ...]
    rationale: str
    cautions: Tuple[str, ...] = ()
    provenance: Dict[str, object] = None  # type: ignore[assignment]
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if self.provenance is None:
            object.__setattr__(self, "provenance", {})

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ConditionRuleSet:
    """Versioned rules plus their deterministic tier-selection policy."""

    definition_id: str
    selection_mode: Literal["first_nonempty_tier"]
    tier_order: Tuple[RuleTier, ...]
    rules: Tuple[ConditionRule, ...]
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


__all__ = [
    "ConditionRule",
    "ConditionRuleSet",
    "PartnerConstraint",
    "RuleRecommendation",
    "RuleScope",
    "RuleSelection",
    "RuleStatus",
    "RuleTier",
]
