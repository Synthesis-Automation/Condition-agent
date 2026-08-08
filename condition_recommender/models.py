"""Versioned records shared by conversion and future recommendation stages."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, Literal, Optional, Tuple


RECOMMENDATION_RECORD_SCHEMA_VERSION = "10.1"
GENERIC_CONVERTER_DEFINITION_VERSION = "generic_conversion.v10.1"
COMPATIBLE_RECOMMENDATION_RECORD_SCHEMA_VERSIONS = frozenset({"10.0", "10.1"})
COMPATIBLE_GENERIC_CONVERTER_DEFINITION_VERSIONS = frozenset(
    {"generic_conversion.v10.0", "generic_conversion.v10.1"}
)
CORE_ELIGIBILITY_DEFINITION_VERSION = "core_eligibility.v1@1.0"
CHEMIST_RANKING_PREFERENCES_SCHEMA_VERSION = "1.0"
GENERIC_RECOMMENDATION_RESULT_SCHEMA_VERSION = "3.4"
REACTION_COMPLETION_PROPOSAL_SCHEMA_VERSION = "1.0"
FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION = (
    "fragment_source_capabilities.v1@1.2"
)


class AdmissionTier(str, Enum):
    VERIFIED = "verified"
    REVIEW = "review"
    REJECTED = "rejected"


class ChemistryStatus(str, Enum):
    UNCLASSIFIED = "unclassified"
    VERIFIED = "verified"
    REVIEW = "review"
    REJECTED = "rejected"


class ConditionStatus(str, Enum):
    UNCLASSIFIED = "unclassified"
    RESOLVED_COMPLETE = "resolved_complete"
    RESOLVED_PARTIAL = "resolved_partial"
    UNRESOLVED_RETAINED = "unresolved_retained"
    UNUSABLE = "unusable"


class ConditionStageStatus(str, Enum):
    """Whether source conditions are assigned to an ordered process stage."""

    UNCLASSIFIED = "unclassified"
    SINGLE_STAGE = "single_stage"
    STRUCTURED_MULTISTAGE = "structured_multistage"
    UNASSIGNED_MULTISTAGE = "unassigned_multistage"


class OutcomeStatus(str, Enum):
    UNCLASSIFIED = "unclassified"
    USABLE = "usable"
    MISSING = "missing"
    INVALID = "invalid"


class IndexEligibility(str, Enum):
    UNCLASSIFIED = "unclassified"
    ELIGIBLE = "eligible"
    REVIEW_ONLY = "review_only"
    INELIGIBLE = "ineligible"


class CoreEligibility(str, Enum):
    """Trust tier for using a reaction core in recommendation."""

    TRUSTED_CORE = "trusted_core"
    REVIEW_CORE = "review_core"
    QUERY_ONLY = "query_only"
    BLOCKED = "blocked"
    UNAVAILABLE = "unavailable"


class PrecedentTier(str, Enum):
    """Authority granted to one indexed condition precedent."""

    TRUSTED = "trusted"
    REVIEW_CORE = "review_core"


class PrecedentIndexScope(str, Enum):
    """Precedent tiers contained by one persisted retrieval index."""

    TRUSTED = "trusted"
    TRUSTED_AND_REVIEW_CORE = "trusted_and_review_core"


@dataclass(frozen=True)
class ChemistRankingPreferences:
    """Versioned chemist-selected priorities for recipe reranking.

    Chemistry admission, reaction compatibility, and hard condition conflicts
    are intentionally outside this contract and cannot be disabled by weights.
    An empty ``weights`` mapping requests the named declarative preset.
    """

    profile_id: str = "default"
    weights: Dict[str, float] = field(default_factory=dict)
    definition_id: str = "chemist_ranking_profiles.v1"
    definition_version: str = "1.0"
    customized: bool = False
    schema_version: str = CHEMIST_RANKING_PREFERENCES_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


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

        return "COND1:" + "|".join(
            (
                token("catalyst", self.catalyst_cas),
                token("reagent", self.reagent_cas),
                token("solvent", self.solvent_cas),
            )
        )


@dataclass(frozen=True)
class ReferenceIdentity:
    """Deterministic publication identity with source-faithful provenance."""

    reference_id: str
    raw_reference: str
    normalized_citation: str
    doi: Optional[str] = None
    patent_number: Optional[str] = None
    publication_year: Optional[int] = None
    resolution_status: str = "missing"
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

@dataclass(frozen=True)
class FragmentSourceSupport:
    """Condition evidence that can supply one product-only fragment."""

    requirement_id: str
    fragment_key: str
    status: str
    component_substance_ids: Tuple[str, ...] = ()
    component_raw_identifiers: Tuple[str, ...] = ()
    capability_ids: Tuple[str, ...] = ()
    evidence: Tuple[str, ...] = ()
    confidence: float = 0.0
    definition_version: str = FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ReactionCompletionOption:
    """One curated way to satisfy a product-fragment source requirement."""

    option_id: str
    option_kind: Literal[
        "compatible_source_class",
        "registered_substance",
        "unresolved",
    ]
    display_name: str
    capability_id: Optional[str] = None
    substance_id: Optional[str] = None
    canonical_name: Optional[str] = None
    schema_version: str = REACTION_COMPLETION_PROPOSAL_SCHEMA_VERSION


@dataclass(frozen=True)
class ReactionCompletionRequirement:
    """A structural source gap and the curated choices that may satisfy it."""

    requirement_id: str
    fragment_key: str
    canonical_fragment_smiles: str
    rooted_fragment_smiles: str
    element_counts: Dict[str, int]
    attachment_element: str
    options: Tuple[ReactionCompletionOption, ...]
    schema_version: str = REACTION_COMPLETION_PROPOSAL_SCHEMA_VERSION


@dataclass(frozen=True)
class ReactionCompletionProposal:
    """System-proposed condition-source completion for an incomplete query."""

    query_reaction_smiles: str
    proposal_id: str
    status: Literal[
        "not_required",
        "confirmation_recommended",
        "no_curated_source_options",
    ]
    requirements: Tuple[ReactionCompletionRequirement, ...] = ()
    warnings: Tuple[str, ...] = ()
    provenance: Literal["system_proposed"] = "system_proposed"
    definition_version: str = FRAGMENT_SOURCE_CAPABILITY_DEFINITION_VERSION
    schema_version: str = REACTION_COMPLETION_PROPOSAL_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ReactionCompletionSelection:
    """A user decision about one system-proposed source requirement."""

    proposal_id: str
    requirement_id: str
    selection_kind: Literal[
        "compatible_source_class",
        "registered_substance",
        "custom_identifier",
        "unresolved",
    ]
    provenance: Literal["user_confirmed", "user_edited"]
    display_name: str
    capability_id: Optional[str] = None
    substance_id: Optional[str] = None
    raw_identifier: Optional[str] = None
    resolved: bool = False
    schema_version: str = REACTION_COMPLETION_PROPOSAL_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.proposal_id or not self.requirement_id:
            raise ValueError("completion selection requires proposal and requirement IDs")
        if (
            self.selection_kind == "compatible_source_class"
            and not self.capability_id
        ):
            raise ValueError("compatible source selection requires a capability ID")
        if self.selection_kind == "registered_substance" and (
            not self.capability_id or not self.substance_id or not self.resolved
        ):
            raise ValueError(
                "registered source selection requires resolved capability and substance IDs"
            )
        if self.selection_kind == "custom_identifier" and not self.raw_identifier:
            raise ValueError("custom source selection requires a raw identifier")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


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
    reaction_label: Optional[Dict[str, Any]]
    yield_pct: Optional[float]
    temperature_c: Optional[float]
    time_h: Optional[float]
    conditions: ConditionIdentity
    chemistry_status: ChemistryStatus = ChemistryStatus.UNCLASSIFIED
    condition_status: ConditionStatus = ConditionStatus.UNCLASSIFIED
    condition_stage_status: ConditionStageStatus = ConditionStageStatus.UNCLASSIFIED
    outcome_status: OutcomeStatus = OutcomeStatus.UNCLASSIFIED
    index_eligibility: IndexEligibility = IndexEligibility.UNCLASSIFIED
    precedent_tier: Optional[PrecedentTier] = None
    core_eligibility: CoreEligibility = CoreEligibility.UNAVAILABLE
    core_eligibility_reasons: Tuple[str, ...] = ()
    core_eligibility_definition_version: str = CORE_ELIGIBILITY_DEFINITION_VERSION
    family_environment: Optional[Dict[str, Any]] = None
    product_connection: Optional[Dict[str, Any]] = None
    spectator_groups: Tuple[Dict[str, Any], ...] = ()
    partial_product_transformation: Optional[Dict[str, Any]] = None
    reaction_completeness: Optional[Dict[str, Any]] = None
    reaction_signature: Optional[Dict[str, Any]] = None
    reaction_core: Optional[Dict[str, Any]] = None
    reaction_evidence_candidates: Tuple[Dict[str, Any], ...] = ()
    reaction_edit_hypotheses: Tuple[Dict[str, Any], ...] = ()
    external_atom_mapping: Optional[Dict[str, Any]] = None
    fallback_descriptor: Optional[Dict[str, Any]] = None
    reaction_observation: Optional[Dict[str, Any]] = None
    reaction_interpretation: Optional[Dict[str, Any]] = None
    fragment_source_support: Tuple[Dict[str, Any], ...] = ()
    transformation_class: Optional[str] = None
    transformation_confidence: float = 0.0
    family_confidence: float = 0.0
    taxonomy_definition_versions: Dict[str, str] = field(default_factory=dict)
    source_dataset: str = ""
    source_path: str = ""
    source_declared_family: str = ""
    reference_id: str = ""
    reference_identity: Optional[Dict[str, Any]] = None
    observation_id: str = ""
    canonical_reaction_id: Optional[str] = None
    canonical_reaction_smiles: Optional[str] = None
    raw_recipe_id: str = ""
    resolved_recipe_core_id: str = ""
    condition_resolution: Dict[str, Any] = field(default_factory=dict)
    resolved_recipe_id: str = ""
    resolved_recipe: Optional[Dict[str, Any]] = None
    synthesis_protocol: Optional[Dict[str, Any]] = None
    reference_condition_series_id: str = ""
    converter_definition_version: str = GENERIC_CONVERTER_DEFINITION_VERSION
    source: Dict[str, Any] = field(default_factory=dict)
    schema_version: str = RECOMMENDATION_RECORD_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GenericConditionRecommendation:
    """One canonical recipe aggregated from admitted structural precedents."""

    rank: int
    recipe_id: str
    recipe_core_id: str
    recipe_variant_ids: Tuple[str, ...]
    resolved_recipe: Dict[str, Any]
    score: float
    similarity_score: float
    compatibility_score: float
    expected_yield_pct: Optional[float]
    support: int
    observation_support: int
    reference_support: int
    condition_series_support: int
    dataset_support: int
    retrieval_level: str
    precedent_reaction_ids: Tuple[str, ...]
    precedent_reference_ids: Tuple[str, ...]
    explanation: Tuple[str, ...]
    score_trace: "RecommendationScoreTrace"
    synthesis_protocol: Dict[str, Any] = field(default_factory=dict)
    default_rank: Optional[int] = None
    default_score: Optional[float] = None
    rank_change: int = 0
    factor_evidence: Dict[str, Any] = field(default_factory=dict)
    precedent_reaction_smiles: Tuple[str, ...] = ()
    precedent_reaction_contexts: Tuple[Dict[str, Any], ...] = ()
    compatibility_evidence: Tuple[str, ...] = ()
    cautions: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RetrievalLevelTrace:
    """One attempted retrieval tier with correlated support made explicit."""

    level: str
    candidate_count: int
    independent_candidate_count: int
    compatible_candidate_count: int
    independent_compatible_candidate_count: int
    excluded_candidate_count: int
    minimum_independent_support: int
    status: str


@dataclass(frozen=True)
class RecommendationScoreTrace:
    """Auditable similarity, ranking, support, and outcome contributions."""

    similarity_components: Dict[str, float]
    similarity_contributions: Dict[str, float]
    ranking_components: Dict[str, Optional[float]]
    ranking_contributions: Dict[str, float]
    applied_ranking_weights: Dict[str, float]
    independent_evidence_count: int
    observed_outcome_count: int
    pool_yield_prior_pct: Optional[float]
    definition_versions: Dict[str, str]
    ranking_profile: str = "default"
    default_ranking_contributions: Dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class GenericRecommendationResult:
    """Type-agnostic recommendation result over generic converted records."""

    query_reaction_smiles: str
    valid: bool
    effective_query_reaction_smiles: Optional[str] = None
    query_signature_id: Optional[str] = None
    query_reaction_core_id: Optional[str] = None
    query_fallback_descriptor_id: Optional[str] = None
    query_edit_hypothesis_ids: Tuple[str, ...] = ()
    external_mapping_status: Optional[str] = None
    external_mapping_provider: Optional[str] = None
    external_mapping_confidence: Optional[float] = None
    recommendation_mode: str = "verified_signature"
    reaction_label: Optional[Dict[str, Any]] = None
    named_family: Optional[str] = None
    transformation_class: Optional[str] = None
    spectator_groups: Tuple[Dict[str, Any], ...] = ()
    reaction_partners: Tuple[Dict[str, Any], ...] = ()
    completion_proposal: Optional[Dict[str, Any]] = None
    completion_selections: Tuple[Dict[str, Any], ...] = ()
    ranking_preferences: Dict[str, Any] = field(default_factory=dict)
    retrieval_definition_version: str = ""
    retrieval_strategy: str = "hybrid"
    retrieval_level: Optional[str] = None
    candidate_count: int = 0
    independent_candidate_count: int = 0
    compatible_candidate_count: int = 0
    independent_compatible_candidate_count: int = 0
    excluded_candidate_count: int = 0
    retrieval_trace: Tuple[RetrievalLevelTrace, ...] = ()
    recommendations: Tuple[GenericConditionRecommendation, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = GENERIC_RECOMMENDATION_RESULT_SCHEMA_VERSION

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
