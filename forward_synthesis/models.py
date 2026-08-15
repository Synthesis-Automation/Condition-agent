"""Immutable public contracts for forward reaction prediction."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Optional, Tuple

from reactive_taxonomy import BidirectionalReactionOperator, OperatorAtomCorrespondence

from .condition_profiles import (
    ForwardConditionProfile,
    ForwardConditionProfileEvidence,
)


FORWARD_LIBRARY_SCHEMA_VERSION = "1.0"
FORWARD_PREDICTION_SCHEMA_VERSION = "1.2"
FORWARD_ROUTE_ASSESSMENT_SCHEMA_VERSION = "1.2"

ForwardDisposition = Literal[
    "clear",
    "competitive",
    "unsupported",
    "structurally_inconsistent",
    "out_of_scope",
]


@dataclass(frozen=True)
class ForwardPrecursorIndex:
    """Conservative necessary-feature index for precursor-side retrieval."""

    component_count_to_operator_ids: Dict[int, Tuple[str, ...]]
    operator_required_atomic_numbers: Dict[str, Tuple[int, ...]]
    definition_id: str = "forward_precursor_index.v1"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardOperatorLibrary:
    """Versioned collection of source-forward-round-tripped operators."""

    operators: Tuple[BidirectionalReactionOperator, ...]
    source_template_count: int
    admitted_operator_count: int
    rejection_counts: Dict[str, int]
    precursor_index: ForwardPrecursorIndex
    definition_id: str = "forward_operator_library.v1"
    source_library_definition_id: str = ""
    schema_version: str = FORWARD_LIBRARY_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != FORWARD_LIBRARY_SCHEMA_VERSION:
            raise ValueError("unsupported forward operator-library schema")
        ids = [operator.forward_operator_id for operator in self.operators]
        if len(ids) != len(set(ids)):
            raise ValueError("forward operator IDs must be unique")
        if self.admitted_operator_count != len(self.operators):
            raise ValueError("admitted operator count contradicts library contents")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardRecipeEvidence:
    """Compatibility of one supplied recipe with one generated reaction."""

    evaluated: bool
    compatible: Optional[bool]
    score: Optional[float]
    hard_conflicts: Tuple[str, ...] = ()
    cautions: Tuple[str, ...] = ()
    definition_id: str = ""
    definition_version: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardProductCandidate:
    """One structurally validated predicted product and its evidence trace."""

    rank: int
    product_smiles: str
    reaction_smiles: str
    mapped_reaction_smiles: str
    pathway_id: str
    operator_id: str
    forward_operator_id: str
    realization_id: str
    template_id: str
    abstraction_level: str
    participating_component_indices: Tuple[int, ...]
    participating_precursor_smiles: str
    assignment: Tuple[int, ...]
    reactant_stoichiometry: Tuple[Tuple[int, int], ...]
    uses_virtual_copies: bool
    atom_correspondence: Tuple[OperatorAtomCorrespondence, ...]
    score: float
    score_components: Dict[str, float]
    structural_score_band: int
    reverse_round_trip: bool
    reaction_signature_id: Optional[str]
    reaction_signature_schema_version: Optional[str]
    operator_edit_agreement: bool
    observed_edit_tokens: Tuple[str, ...]
    independent_reference_support: int
    observation_support: int
    precedent_reaction_ids: Tuple[str, ...]
    precedent_reference_ids: Tuple[str, ...]
    recipe_evidence: ForwardRecipeEvidence
    condition_profile_evidence: ForwardConditionProfileEvidence
    alternative_pathway_ids: Tuple[str, ...] = ()
    alternative_operator_ids: Tuple[str, ...] = ()
    alternative_template_ids: Tuple[str, ...] = ()
    named_annotations: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    schema_version: str = FORWARD_PREDICTION_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardSearchDiagnostics:
    """Counts for every precursor retrieval and graph-validation stage."""

    library_operator_count: int = 0
    indexed_operator_count: int = 0
    applied_operator_count: int = 0
    generated_outcome_count: int = 0
    reverse_round_trip_failure_count: int = 0
    invalid_reaction_count: int = 0
    missing_signature_count: int = 0
    operator_edit_mismatch_count: int = 0
    recipe_conflict_count: int = 0
    valid_pathway_count: int = 0
    self_reaction_pathway_count: int = 0
    unique_product_count: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardCompetitionGroup:
    """Products competing through one operator or one concrete site assignment."""

    competition_level: Literal["operator", "site", "product"]
    group_key: str
    candidate_ranks: Tuple[int, ...]
    product_smiles: Tuple[str, ...]
    operator_ids: Tuple[str, ...]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ForwardPredictionResult:
    """Blind product prediction from starting materials and optional conditions."""

    query_starting_materials: str
    canonical_starting_materials: str
    conditions_supplied: bool
    condition_profile_supplied: bool
    condition_profile: ForwardConditionProfile
    self_reactions_considered: bool
    valid: bool
    status: str
    candidates: Tuple[ForwardProductCandidate, ...]
    competition_groups: Tuple[ForwardCompetitionGroup, ...]
    diagnostics: ForwardSearchDiagnostics
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    ranking_definition_id: str = "forward_ranking.v3"
    schema_version: str = FORWARD_PREDICTION_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RouteStepForwardAssessment:
    """Targeted replay plus target-blind competition audit for a planned step."""

    starting_materials: str
    intended_product: str
    intended_match: Literal[
        "exact", "stereo_relaxed", "connectivity_only", "absent", "invalid"
    ]
    targeted_replay_status: str
    intended_product_rank: Optional[int]
    best_competitor_product: Optional[str]
    score_margin: Optional[float]
    disposition: ForwardDisposition
    blind_prediction: ForwardPredictionResult
    operator_hint: Optional[str] = None
    warnings: Tuple[str, ...] = ()
    schema_version: str = FORWARD_ROUTE_ASSESSMENT_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


__all__ = [
    "FORWARD_LIBRARY_SCHEMA_VERSION",
    "FORWARD_PREDICTION_SCHEMA_VERSION",
    "FORWARD_ROUTE_ASSESSMENT_SCHEMA_VERSION",
    "ForwardCompetitionGroup",
    "ForwardDisposition",
    "ForwardOperatorLibrary",
    "ForwardPrecursorIndex",
    "ForwardPredictionResult",
    "ForwardProductCandidate",
    "ForwardRecipeEvidence",
    "ForwardSearchDiagnostics",
    "RouteStepForwardAssessment",
]
