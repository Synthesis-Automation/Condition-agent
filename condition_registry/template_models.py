"""Typed contracts for versioned expert condition-recipe templates."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Literal, Optional, Tuple


TemplateStatus = Literal["draft", "active", "retired"]


@dataclass(frozen=True)
class RecipeTemplateOption:
    """One registry-backed choice for a recipe slot."""

    substance_id: str
    preference: int = 0
    notes: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RecipeTemplateSlot:
    """Role-constrained alternatives without implicit Cartesian expansion."""

    slot_id: str
    role_id: str
    required: bool
    alternatives: Tuple[RecipeTemplateOption, ...]
    selection_policy: Literal["present_alternatives", "select_one"] = (
        "present_alternatives"
    )
    notes: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RecipeTemplateVariant:
    """One explicit slot selection; variants are never generated implicitly."""

    variant_id: str
    status: TemplateStatus
    selections: Tuple["RecipeTemplateSelection", ...]
    priority: int = 0
    notes: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RecipeTemplateSelection:
    """Registry identity and quantity selected for one recipe slot."""

    slot_id: str
    substance_id: str
    amount: Optional[float] = None
    amount_unit: Optional[str] = None

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RecipeTemplatePartnerAmount:
    """Stoichiometric range for one structurally assigned reaction partner."""

    role: str
    minimum: float
    maximum: float
    unit: Literal["equivalent"] = "equivalent"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ConditionRecipeTemplate:
    """Versioned expert recipe whose incomplete details remain explicit."""

    template_id: str
    status: TemplateStatus
    slots: Tuple[RecipeTemplateSlot, ...]
    variants: Tuple[RecipeTemplateVariant, ...] = ()
    temperature_c: Optional[float] = None
    time_h: Optional[float] = None
    concentration_m: Optional[float] = None
    pressure_bar: Optional[float] = None
    atmosphere: Optional[str] = None
    partner_amounts: Tuple[RecipeTemplatePartnerAmount, ...] = ()
    forbidden_substance_ids: Tuple[str, ...] = ()
    notes: Tuple[str, ...] = ()
    provenance: Dict[str, object] = None  # type: ignore[assignment]
    schema_version: str = "1.2"

    def __post_init__(self) -> None:
        if self.provenance is None:
            object.__setattr__(self, "provenance", {})

    @property
    def identity_complete(self) -> bool:
        """Return whether every required choice is registry-backed."""
        return all(
            slot.alternatives
            and all(option.substance_id for option in slot.alternatives)
            for slot in self.slots
            if slot.required
        )

    def to_dict(self) -> Dict[str, object]:
        payload = asdict(self)
        payload["identity_complete"] = self.identity_complete
        return payload


@dataclass(frozen=True)
class ConditionRecipeTemplateSet:
    """One immutable collection of condition-recipe templates."""

    definition_id: str
    templates: Tuple[ConditionRecipeTemplate, ...]
    schema_version: str = "1.2"

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


__all__ = [
    "ConditionRecipeTemplate",
    "ConditionRecipeTemplateSet",
    "RecipeTemplateOption",
    "RecipeTemplatePartnerAmount",
    "RecipeTemplateSelection",
    "RecipeTemplateSlot",
    "RecipeTemplateVariant",
    "TemplateStatus",
]
