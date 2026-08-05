"""Typed public contracts for condition substance resolution."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Optional, Tuple


CONDITION_IDENTIFIER_TYPES = (
    "canonical_name",
    "common_name",
    "systematic_name",
    "abbreviation",
    "trade_name",
    "legacy_name",
    "cas",
    "inchi_key",
    "database_id",
)

CONDITION_NAME_IDENTIFIER_TYPES = (
    "canonical_name",
    "common_name",
    "systematic_name",
    "abbreviation",
    "trade_name",
    "legacy_name",
)

CONDITION_RECIPE_COMPONENT_BUCKETS = (
    "catalysts",
    "ligands",
    "bases",
    "acids",
    "condensation_agents",
    "oxidants",
    "reductants",
    "additives",
    "solvents",
    "other_components",
)


@dataclass(frozen=True)
class ConditionComponentInput:
    """Raw condition component awaiting registry identity resolution."""

    raw_identifier: str
    source_field: str
    identifier_type: str = "auto"
    source_role_hint: Optional[str] = None
    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        supported = {
            "auto",
            "name",
            "substance_id",
            *CONDITION_IDENTIFIER_TYPES,
        }
        if self.identifier_type not in supported:
            raise ValueError(f"Unsupported condition identifier type: {self.identifier_type}")


@dataclass(frozen=True)
class ConditionProcessStage:
    """One explicit temperature/time segment in a condition procedure."""

    stage_index: int
    temperature_c: Optional[float] = None
    time_h: Optional[float] = None
    atmosphere: Optional[str] = None
    provenance: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class RoleAssignment:
    role_id: str
    family_id: Optional[str]
    tag: Optional[str] = None
    confidence: float = 1.0
    evidence: str = "curated_registry"


@dataclass(frozen=True)
class SubstanceIdentifier:
    """One typed, provenance-bearing identifier for a substance."""

    identifier_id: str
    substance_id: str
    identifier_type: str
    value: str
    language: Optional[str] = None
    is_preferred: bool = False
    source: Optional[str] = None
    confidence: float = 1.0
    status: str = "active"
    normalization_profile: Optional[str] = None
    allow_ambiguous: bool = False

    def __post_init__(self) -> None:
        if self.identifier_type not in CONDITION_IDENTIFIER_TYPES:
            raise ValueError(
                f"Unsupported substance identifier type: {self.identifier_type}"
            )
        if not self.identifier_id.strip():
            raise ValueError("Substance identifier ID must not be empty")
        if not self.substance_id.strip():
            raise ValueError("Substance identifier must reference a substance")
        if not self.value.strip():
            raise ValueError("Substance identifier value must not be empty")
        if self.status not in {"active", "deprecated"}:
            raise ValueError(f"Unsupported substance identifier status: {self.status}")
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("Substance identifier confidence must be between 0 and 1")


@dataclass(frozen=True)
class Substance:
    substance_id: str
    canonical_name: str
    cas: Optional[str]
    smiles: Optional[str]
    aliases: Tuple[str, ...] = ()
    roles: Tuple[RoleAssignment, ...] = ()
    properties: Dict[str, Any] = None  # type: ignore[assignment]
    identifiers: Tuple[SubstanceIdentifier, ...] = ()

    def __post_init__(self) -> None:
        if self.properties is None:
            object.__setattr__(self, "properties", {})
        identifiers = list(self.identifiers)
        identifier_values = {
            identifier.value for identifier in identifiers
        }
        for index, alias in enumerate(self.aliases, start=1):
            if alias in identifier_values:
                continue
            identifiers.append(
                SubstanceIdentifier(
                    identifier_id=(
                        f"compat:{self.substance_id}:legacy_alias:{index}"
                    ),
                    substance_id=self.substance_id,
                    identifier_type="legacy_name",
                    value=alias,
                    source="Substance.aliases compatibility field",
                    normalization_profile="chemical_name_v1",
                )
            )
        object.__setattr__(self, "identifiers", tuple(identifiers))
        aliases = tuple(
            dict.fromkeys(
                identifier.value
                for identifier in identifiers
                if identifier.status == "active"
                and identifier.identifier_type in CONDITION_NAME_IDENTIFIER_TYPES
                and identifier.identifier_type != "canonical_name"
            )
        )
        object.__setattr__(self, "aliases", aliases)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ResolutionResult:
    query: str
    status: str
    substance: Optional[Substance] = None
    match_kind: Optional[str] = None
    candidates: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    matched_identifier: Optional[SubstanceIdentifier] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ContextualRoleAssignment:
    """One context-ranked role interpretation with explicit evidence."""

    role_id: str
    family_id: Optional[str]
    confidence: float
    evidence: Tuple[str, ...]


@dataclass(frozen=True)
class ResolvedConditionComponent:
    """Condition identity and contextual roles retained with provenance."""

    raw_identifier: str
    source_field: str
    identity_status: str
    substance_id: Optional[str]
    canonical_name: Optional[str]
    roles: Tuple[ContextualRoleAssignment, ...]
    primary_role: str
    primary_role_confidence: float
    amount: Optional[float] = None
    amount_unit: Optional[str] = None
    source_role_hint: Optional[str] = None
    warnings: Tuple[str, ...] = ()
    provenance: Dict[str, Any] = field(default_factory=dict)
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ResolvedConditionRecipe:
    """Canonical role-aware recipe used to group equivalent precedents."""

    recipe_id: str
    recipe_core_id: str
    catalysts: Tuple[ResolvedConditionComponent, ...] = ()
    ligands: Tuple[ResolvedConditionComponent, ...] = ()
    bases: Tuple[ResolvedConditionComponent, ...] = ()
    acids: Tuple[ResolvedConditionComponent, ...] = ()
    condensation_agents: Tuple[ResolvedConditionComponent, ...] = ()
    oxidants: Tuple[ResolvedConditionComponent, ...] = ()
    reductants: Tuple[ResolvedConditionComponent, ...] = ()
    additives: Tuple[ResolvedConditionComponent, ...] = ()
    solvents: Tuple[ResolvedConditionComponent, ...] = ()
    other_components: Tuple[ResolvedConditionComponent, ...] = ()
    temperature_c: Optional[float] = None
    time_h: Optional[float] = None
    concentration_m: Optional[float] = None
    atmosphere: Optional[str] = None
    stages: Tuple[ConditionProcessStage, ...] = ()
    declared_absences: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    definition_versions: Dict[str, str] = field(default_factory=dict)
    schema_version: str = "1.2"

    @property
    def components(self) -> Tuple[ResolvedConditionComponent, ...]:
        return tuple(
            component
            for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS
            for component in getattr(self, bucket)
        )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
