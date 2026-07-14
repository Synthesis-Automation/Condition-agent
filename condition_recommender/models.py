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
    reaction_signature: Optional[Dict[str, Any]] = None
    transformation_class: Optional[str] = None
    transformation_confidence: float = 0.0
    family_confidence: float = 0.0
    taxonomy_definition_versions: Dict[str, str] = field(default_factory=dict)
    source: Dict[str, Any] = field(default_factory=dict)
    schema_version: str = "1.2"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ConditionRecommendation:
    rank: int
    conditions: ConditionIdentity
    recipe_id: str
    score: float
    similarity_score: float
    expected_yield_pct: float
    support: int
    retrieval_level: str
    precedent_reaction_ids: Tuple[str, ...]
    explanation: Tuple[str, ...]
    cautions: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RecommendationResult:
    query_reaction_smiles: str
    valid: bool
    query_label: Optional[str] = None
    product_connection_label: Optional[str] = None
    retrieval_level: Optional[str] = None
    candidate_count: int = 0
    recommendations: Tuple[ConditionRecommendation, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class LabelConditionRecommendation:
    """One condition recipe retrieved from weak reaction-label precedents."""

    rank: int
    recipe_id: str
    score: float
    label_similarity: float
    signature_similarity: float
    qualifier_similarity: float
    expected_yield_pct: float
    mean_z_score: float
    support: int
    source_reaction_types: Tuple[str, ...]
    source_row_numbers: Tuple[int, ...]
    conditions: Dict[str, str]
    explanation: Tuple[str, ...]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class LabelRecommendationResult:
    """Top condition recipes for a structurally verified reaction query."""

    query_reaction_smiles: str
    valid: bool
    query_label: Optional[str] = None
    grammar_id: Optional[str] = None
    query_signatures: Tuple[str, ...] = ()
    candidate_count: int = 0
    recipe_count: int = 0
    recommendations: Tuple[LabelConditionRecommendation, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
