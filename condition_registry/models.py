"""Typed public contracts for condition substance resolution."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Optional, Tuple


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
        if self.identifier_type not in {"auto", "cas", "name", "substance_id"}:
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
class Substance:
    substance_id: str
    canonical_name: str
    cas: Optional[str]
    smiles: Optional[str]
    aliases: Tuple[str, ...] = ()
    roles: Tuple[RoleAssignment, ...] = ()
    properties: Dict[str, Any] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.properties is None:
            object.__setattr__(self, "properties", {})

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
        return (
            self.catalysts
            + self.ligands
            + self.bases
            + self.acids
            + self.condensation_agents
            + self.oxidants
            + self.reductants
            + self.additives
            + self.solvents
            + self.other_components
        )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
