"""Versioned records shared by conversion and future recommendation stages."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, Optional, Tuple


class AdmissionTier(str, Enum):
    VERIFIED = "verified"
    REVIEW = "review"
    REJECTED = "rejected"


@dataclass(frozen=True)
class ConditionIdentity:
    """Source-faithful condition identities without inferred chemical roles."""

    catalyst_cas: Tuple[str, ...] = ()
    reagent_cas: Tuple[str, ...] = ()
    solvent_cas: Tuple[str, ...] = ()

    @property
    def complete(self) -> bool:
        return bool(self.catalyst_cas and self.reagent_cas and self.solvent_cas)

    @property
    def recipe_id(self) -> str:
        def token(name: str, values: Tuple[str, ...]) -> str:
            return f"{name}={'+'.join(values) if values else 'unknown'}"
        return "COND1:" + "|".join((
            token("catalyst", self.catalyst_cas),
            token("reagent", self.reagent_cas),
            token("solvent", self.solvent_cas),
        ))


@dataclass(frozen=True)
class RecommendationRecord:
    """One converted reaction observation with admission provenance."""

    reaction_id: str
    source_row_number: int
    reaction_smiles: str
    admission_tier: AdmissionTier
    admission_reasons: Tuple[str, ...]
    evidence_quality: str
    named_family: Optional[str]
    reaction_label: Optional[str]
    reaction_label_status: str
    yield_pct: Optional[float]
    temperature_c: Optional[float]
    time_h: Optional[float]
    conditions: ConditionIdentity
    family_environment: Optional[Dict[str, Any]] = None
    product_connection: Optional[Dict[str, Any]] = None
    spectator_groups: Tuple[Dict[str, Any], ...] = ()
    source: Dict[str, Any] = field(default_factory=dict)
    schema_version: str = "1.1"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
