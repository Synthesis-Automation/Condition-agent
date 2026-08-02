"""Typed, source-faithful contracts for chemistry-free dataset ingestion."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Optional, Tuple


INTERMEDIATE_OBSERVATION_SCHEMA_VERSION = "source_observation.v1"


@dataclass(frozen=True)
class SourceIdentifier:
    """One source-supplied identifier for a condition component."""

    identifier_type: str
    value: str
    source_field: str


@dataclass(frozen=True)
class ConditionAmountInput:
    """A normalized amount paired with its source representation."""

    value: Optional[float]
    unit: str
    raw_value: str


@dataclass(frozen=True)
class ConditionComponentGroup:
    """A source slot that may describe a multi-component formulation."""

    group_key: str
    source_slot: str
    display_name: str = ""
    amount: Optional[ConditionAmountInput] = None
    provenance: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ConditionComponentClaim:
    """Source evidence about one condition component before registry resolution."""

    component_key: str
    source_slot: str
    source_role_hint: Optional[str]
    identifiers: Tuple[SourceIdentifier, ...]
    amount: Optional[ConditionAmountInput] = None
    group_key: Optional[str] = None
    introduced_in_stage: Optional[int] = 1
    provenance: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ConditionStageInput:
    """One source-observed process stage without chemical interpretation."""

    stage_index: int
    component_keys: Tuple[str, ...] = ()
    temperature_c: Optional[float] = None
    time_h: Optional[float] = None
    atmosphere: Optional[str] = None
    pressure_bar: Optional[float] = None
    concentration_m: Optional[float] = None
    source_text: str = ""
    provenance: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ConditionInput:
    """All source-observed material and operating-condition claims."""

    components: Tuple[ConditionComponentClaim, ...] = ()
    component_groups: Tuple[ConditionComponentGroup, ...] = ()
    stages: Tuple[ConditionStageInput, ...] = ()
    declared_stage_count: Optional[int] = None
    procedure_text: str = ""
    declared_absences: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class ReactionEvidenceInput:
    """Either source molecular structure or explicitly unverified labels."""

    evidence_kind: str
    reaction_smiles: Optional[str] = None
    supplied_mapping_status: str = "not_applicable"
    source_reaction_type: str = ""
    source_labels: Dict[str, Any] = field(default_factory=dict)
    structure_available: bool = False


@dataclass(frozen=True)
class OutcomeInput:
    """One reported outcome with its measurement basis retained."""

    outcome_type: str
    value: Optional[float]
    unit: str
    raw_value: str
    source_field: str
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class SourceProvenance:
    """Stable source location and experimental grouping metadata."""

    corpus_id: str
    release_id: str
    adapter_id: str
    adapter_version: str
    source_file: str
    source_file_sha256: str
    source_row_number: int
    source_record_id: str
    source_groups: Dict[str, str] = field(default_factory=dict)
    reference: str = ""


@dataclass(frozen=True)
class CanonicalSourceObservation:
    """Uniform chemistry-free observation produced from one source row."""

    observation_id: str
    observation_kind: str
    source: SourceProvenance
    reaction: ReactionEvidenceInput
    conditions: ConditionInput
    outcomes: Tuple[OutcomeInput, ...] = ()
    ingestion_status: str = "accepted"
    warnings: Tuple[str, ...] = ()
    raw_fields: Dict[str, Any] = field(default_factory=dict)
    schema_version: str = INTERMEDIATE_OBSERVATION_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable nested representation."""
        return asdict(self)


__all__ = [
    "INTERMEDIATE_OBSERVATION_SCHEMA_VERSION",
    "CanonicalSourceObservation",
    "ConditionAmountInput",
    "ConditionComponentClaim",
    "ConditionComponentGroup",
    "ConditionInput",
    "ConditionStageInput",
    "OutcomeInput",
    "ReactionEvidenceInput",
    "SourceIdentifier",
    "SourceProvenance",
]
