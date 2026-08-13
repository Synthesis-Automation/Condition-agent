"""Contracts for structurally diverse core-derived retrosynthesis templates."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Optional, Tuple

from cas_tools import PrecursorRealismAssessment

from .models import CenterReactivityContext, TemplateContext
from .selectivity_poc import FunctionalGroupCompetitionWarning


GenericLevel = Literal["RDCHIRAL", "L0", "L1", "L2"]


@dataclass(frozen=True)
class GenericTemplatePrecedent:
    """One source-round-tripped precedent for a generic graph operator."""

    reaction_id: str
    reference_id: str
    product_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str
    context: TemplateContext


@dataclass(frozen=True)
class GenericCoreTemplate:
    """One executable template classified only by structural transformation."""

    template_id: str
    operator_id: str
    transformation_kind: Optional[str]
    abstraction_level: GenericLevel
    compiler_engine: Literal["rdchiral", "reaction_core"]
    reaction_smarts: str
    product_smarts: str
    precursor_smarts: str
    edit_tokens: Tuple[str, ...]
    handle_signature: str
    stereo_policy: Literal["exact", "relaxed"]
    observation_support: int
    independent_reference_support: int
    precedents: Tuple[GenericTemplatePrecedent, ...]
    realization_id: str = ""
    operator_signature: str = ""
    synthon_signature: str = ""
    named_annotations: Tuple[str, ...] = ()


@dataclass(frozen=True)
class GenericGraphOperator:
    """One handle-independent graph transformation mined from observations."""

    operator_id: str
    operator_signature: str
    edit_tokens: Tuple[str, ...]
    realization_ids: Tuple[str, ...]
    abstraction_levels: Tuple[str, ...]
    observation_support: int
    independent_reference_support: int
    named_annotations: Tuple[str, ...] = ()


@dataclass(frozen=True)
class GenericHandleCompletionGroup:
    """Data-derived precursor completions sharing an operator and handle class."""

    completion_group_id: str
    operator_id: str
    completion_signature: str
    synthon_signatures: Tuple[str, ...]
    realization_ids: Tuple[str, ...]
    template_ids: Tuple[str, ...]
    handle_signatures: Tuple[str, ...]
    observation_support: int
    independent_reference_support: int


@dataclass(frozen=True)
class GenericRetrievalIndex:
    """Serializable necessary-feature index for product-side template lookup."""

    token_to_template_ids: Dict[str, Tuple[str, ...]]
    template_required_tokens: Dict[str, Tuple[str, ...]]
    fallback_template_ids: Tuple[str, ...] = ()
    definition_id: str = "generic_product_retrieval_index.v1"


@dataclass(frozen=True)
class GenericTemplateLibrary:
    """Serializable generic template collection."""

    templates: Tuple[GenericCoreTemplate, ...]
    source_row_count: int
    accepted_observation_count: int
    rejection_counts: Dict[str, int]
    definition: Dict[str, Any]
    operators: Tuple[GenericGraphOperator, ...] = ()
    completion_groups: Tuple[GenericHandleCompletionGroup, ...] = ()
    retrieval_index: Optional[GenericRetrievalIndex] = None
    schema_version: str = "3.0"

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible representation."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> "GenericTemplateLibrary":
        """Load a generic library and its nested contexts."""

        if value.get("schema_version") not in {"1.0", "2.0", "3.0"}:
            raise ValueError("unsupported generic template library schema")
        templates = []
        for raw_template in value.get("templates") or ():
            precedents = []
            for raw_precedent in raw_template.get("precedents") or ():
                raw_context = raw_precedent.get("context") or {}
                context = TemplateContext(
                    spectator_group_ids=tuple(
                        raw_context.get("spectator_group_ids") or ()
                    ),
                    substituent_feature_tokens=tuple(
                        raw_context.get("substituent_feature_tokens") or ()
                    ),
                    center_profiles=tuple(
                        CenterReactivityContext(**profile)
                        for profile in raw_context.get("center_profiles") or ()
                    ),
                )
                precedents.append(
                    GenericTemplatePrecedent(
                        **{
                            key: item
                            for key, item in raw_precedent.items()
                            if key != "context"
                        },
                        context=context,
                    )
                )
            templates.append(
                GenericCoreTemplate(
                    **{
                        key: item
                        for key, item in raw_template.items()
                        if key
                        not in {
                            "edit_tokens",
                            "precedents",
                            "realization_id",
                            "operator_signature",
                            "synthon_signature",
                            "named_annotations",
                        }
                    },
                    edit_tokens=tuple(raw_template.get("edit_tokens") or ()),
                    precedents=tuple(precedents),
                    realization_id=str(raw_template.get("realization_id") or ""),
                    operator_signature=str(
                        raw_template.get("operator_signature") or ""
                    ),
                    synthon_signature=str(
                        raw_template.get("synthon_signature") or ""
                    ),
                    named_annotations=tuple(
                        raw_template.get("named_annotations") or ()
                    ),
                )
            )
        operators = tuple(
            GenericGraphOperator(
                **{
                    key: item
                    for key, item in raw.items()
                    if key
                    not in {
                        "edit_tokens",
                        "realization_ids",
                        "abstraction_levels",
                        "named_annotations",
                    }
                },
                edit_tokens=tuple(raw.get("edit_tokens") or ()),
                realization_ids=tuple(raw.get("realization_ids") or ()),
                abstraction_levels=tuple(raw.get("abstraction_levels") or ()),
                named_annotations=tuple(raw.get("named_annotations") or ()),
            )
            for raw in value.get("operators") or ()
        )
        completion_groups = tuple(
            GenericHandleCompletionGroup(
                **{
                    key: item
                    for key, item in raw.items()
                    if key
                    not in {
                        "realization_ids",
                        "template_ids",
                        "handle_signatures",
                        "completion_signature",
                        "synthon_signature",
                        "synthon_signatures",
                    }
                },
                completion_signature=str(
                    raw.get("completion_signature")
                    or ";".join(raw.get("handle_signatures") or ())
                    or "legacy_completion"
                ),
                synthon_signatures=tuple(
                    raw.get("synthon_signatures")
                    or (
                        (raw.get("synthon_signature"),)
                        if raw.get("synthon_signature")
                        else ()
                    )
                ),
                realization_ids=tuple(raw.get("realization_ids") or ()),
                template_ids=tuple(raw.get("template_ids") or ()),
                handle_signatures=tuple(raw.get("handle_signatures") or ()),
            )
            for raw in value.get("completion_groups") or ()
        )
        raw_index = value.get("retrieval_index")
        retrieval_index = None
        if isinstance(raw_index, dict):
            retrieval_index = GenericRetrievalIndex(
                token_to_template_ids={
                    str(token): tuple(template_ids or ())
                    for token, template_ids in (
                        raw_index.get("token_to_template_ids") or {}
                    ).items()
                },
                template_required_tokens={
                    str(template_id): tuple(tokens or ())
                    for template_id, tokens in (
                        raw_index.get("template_required_tokens") or {}
                    ).items()
                },
                fallback_template_ids=tuple(
                    raw_index.get("fallback_template_ids") or ()
                ),
                definition_id=str(
                    raw_index.get("definition_id")
                    or "generic_product_retrieval_index.v1"
                ),
            )
        return cls(
            templates=tuple(templates),
            source_row_count=int(value.get("source_row_count") or 0),
            accepted_observation_count=int(
                value.get("accepted_observation_count") or 0
            ),
            rejection_counts={
                str(key): int(count)
                for key, count in (value.get("rejection_counts") or {}).items()
            },
            definition=dict(value.get("definition") or {}),
            operators=operators,
            completion_groups=completion_groups,
            retrieval_index=retrieval_index,
            schema_version="3.0",
        )


@dataclass(frozen=True)
class GenericSearchDiagnostics:
    """Counts for each retrieval, generation, and validation stage."""

    library_template_count: int = 0
    indexed_template_count: int = 0
    metadata_filtered_template_count: int = 0
    product_query_match_count: int = 0
    applied_template_count: int = 0
    generated_precursor_count: int = 0
    validation_attempt_count: int = 0
    valid_candidate_count: int = 0
    invalid_forward_count: int = 0
    unresolved_identity_count: int = 0
    operator_mismatch_count: int = 0

    def to_dict(self) -> Dict[str, int]:
        """Return stage counters as a JSON-compatible mapping."""

        return asdict(self)


@dataclass(frozen=True)
class OperatorLadderDiagnostics:
    """Aggregated counters for a specificity-ladder expansion."""

    levels_attempted: Tuple[str, ...]
    level_diagnostics: Tuple[Tuple[str, GenericSearchDiagnostics], ...]

    @property
    def proposed_action_count(self) -> int:
        """Return the number of cheaply generated precursor actions."""

        return sum(
            diagnostics.generated_precursor_count
            for _, diagnostics in self.level_diagnostics
        )

    @property
    def validation_attempt_count(self) -> int:
        """Return the number of actions subjected to hard graph validation."""

        return sum(
            diagnostics.validation_attempt_count
            for _, diagnostics in self.level_diagnostics
        )

    @property
    def valid_action_count(self) -> int:
        """Return the number of unique validated actions before ladder merging."""

        return sum(
            diagnostics.valid_candidate_count
            for _, diagnostics in self.level_diagnostics
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible profiling record."""

        return {
            "levels_attempted": list(self.levels_attempted),
            "level_diagnostics": {
                level: diagnostics.to_dict()
                for level, diagnostics in self.level_diagnostics
            },
            "proposed_action_count": self.proposed_action_count,
            "validation_attempt_count": self.validation_attempt_count,
            "valid_action_count": self.valid_action_count,
        }


@dataclass(frozen=True)
class GenericDisconnectionCandidate:
    """One structurally validated generic retrosynthesis proposal."""

    target_smiles: str
    precursor_smiles: str
    proposed_reaction_smiles: str
    transformation_kind: Optional[str]
    abstraction_level: GenericLevel
    compiler_engine: str
    template_id: str
    score: float
    context_similarity: float
    product_similarity: float
    precursor_similarity: float
    template_specificity: float
    independent_reference_support: int
    forward_validation_status: str
    center_transition_key: str
    disconnection_site_key: str
    precedent_reaction_ids: Tuple[str, ...]
    operator_id: str = ""
    realization_id: str = ""
    operator_signature: str = ""
    synthon_signature: str = ""
    pre_diversity_rank: int = 0
    diversity_rank: int = 0
    diversity_group_key: Tuple[str, ...] = ()
    structural_score_band: int = 0
    ranking_policy_definition_id: str = ""
    condition_query_reaction_smiles: str = ""
    selectivity_warnings: Tuple[FunctionalGroupCompetitionWarning, ...] = ()
    precursor_realism_score: Optional[float] = None
    precursor_realism_assessments: Tuple[PrecursorRealismAssessment, ...] = ()
    pre_realism_rank: int = 0
    precursor_realism_rank: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-compatible candidate."""

        return asdict(self)
__all__ = [
    "GenericCoreTemplate",
    "GenericDisconnectionCandidate",
    "GenericGraphOperator",
    "GenericRetrievalIndex",
    "GenericSearchDiagnostics",
    "OperatorLadderDiagnostics",
    "GenericHandleCompletionGroup",
    "GenericTemplateLibrary",
    "GenericTemplatePrecedent",
]
