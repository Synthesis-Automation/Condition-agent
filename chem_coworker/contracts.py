"""Public request and response contracts for the condition coworker."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Literal, Optional, Tuple

from condition_recommender import GenericRecommendationResult
from core_retrosynthesis import (
    MultistepRetrosynthesisResult,
    RetrosynthesisConditionEvidence,
    StrategyProposal,
)


@dataclass(frozen=True)
class CompletionChoice:
    """One user-confirmed source choice from a completion proposal."""

    requirement_id: str
    option_id: Optional[str] = None
    custom_identifier: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.requirement_id.strip():
            raise ValueError("completion requirement_id must not be empty")
        if bool(self.option_id) == bool(self.custom_identifier):
            raise ValueError(
                "completion choice requires exactly one option_id or custom_identifier"
            )


ReviewMode = Literal["off", "auto", "always"]
ReviewVerdict = Literal["keep", "downrank", "flag", "needs_information"]


@dataclass(frozen=True)
class ConditionReviewSettings:
    """Bounded settings for the optional post-recommendation LLM review."""

    mode: ReviewMode = "off"
    provider: str = "openai"
    model: str = "gpt-5.6-terra"
    reasoning_effort: str = "medium"
    max_candidates: int = 10
    max_output_tokens: int = 8_000
    apply_order: bool = True

    def __post_init__(self) -> None:
        if self.mode not in {"off", "auto", "always"}:
            raise ValueError("review mode must be off, auto, or always")
        if self.provider not in {"openai", "aliyun"}:
            raise ValueError("review provider must be openai or aliyun")
        if not self.model.strip():
            raise ValueError("review model must not be empty")
        if self.reasoning_effort not in {
            "none",
            "low",
            "medium",
            "high",
            "xhigh",
            "max",
        }:
            raise ValueError("unsupported review reasoning effort")
        if self.max_candidates < 1 or self.max_candidates > 10:
            raise ValueError("review max_candidates must be between 1 and 10")
        if self.max_output_tokens < 256 or self.max_output_tokens > 20_000:
            raise ValueError("review max_output_tokens must be between 256 and 20000")


@dataclass(frozen=True)
class ConditionCandidateReview:
    """One evidence-linked review verdict for an existing recommendation."""

    recipe_id: str
    original_rank: int
    verdict: ReviewVerdict
    issue_codes: Tuple[str, ...]
    evidence_ids: Tuple[str, ...]
    rationale: str
    confidence: float


@dataclass(frozen=True)
class ConditionGroupReview:
    """One LLM-proposed condition strategy with a deterministic representative."""

    group_id: str
    representative_recipe_id: str
    member_recipe_ids: Tuple[str, ...]
    grouping_basis: Tuple[str, ...]
    evidence_ids: Tuple[str, ...]
    rationale: str


@dataclass(frozen=True)
class ConditionReview:
    """Auditable outcome of the optional LLM review stage."""

    status: Literal["skipped", "completed", "failed"]
    provider: str
    model: str
    trigger_reasons: Tuple[str, ...] = ()
    summary: str = ""
    candidates: Tuple[ConditionCandidateReview, ...] = ()
    groups: Tuple[ConditionGroupReview, ...] = ()
    questions: Tuple[str, ...] = ()
    presentation_recipe_ids: Tuple[str, ...] = ()
    warning: Optional[str] = None
    response_id: Optional[str] = None
    input_tokens: int = 0
    output_tokens: int = 0
    provider_attempts: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serialize review evidence and model metadata."""

        return asdict(self)


@dataclass(frozen=True)
class ConditionRequest:
    """One structure-first condition recommendation request."""

    reaction_smiles: str
    top_k: int = 5
    minimum_pool_size: Optional[int] = None
    unrestricted_fallback: bool = False
    ranking_profile: str = "default"
    ranking_weights: Dict[str, float] = field(default_factory=dict)
    completion_choices: Tuple[CompletionChoice, ...] = ()
    preferred_reaction_ids: Tuple[str, ...] = ()
    review: ConditionReviewSettings = field(default_factory=ConditionReviewSettings)

    def __post_init__(self) -> None:
        if not self.reaction_smiles.strip():
            raise ValueError("reaction_smiles must not be empty")
        if self.top_k < 1 or self.top_k > 50:
            raise ValueError("top_k must be between 1 and 50")
        if self.minimum_pool_size is not None and self.minimum_pool_size < 1:
            raise ValueError("minimum_pool_size must be positive")
        if not self.ranking_profile.strip():
            raise ValueError("ranking_profile must not be empty")
        if any(value < 0.0 for value in self.ranking_weights.values()):
            raise ValueError("ranking weights must be non-negative")
        if self.ranking_weights and sum(self.ranking_weights.values()) <= 0.0:
            raise ValueError("at least one ranking weight must be positive")


@dataclass(frozen=True)
class ConditionResponse:
    """A lossless domain result plus a deterministic human-readable summary."""

    request: ConditionRequest
    result: GenericRecommendationResult
    answer: str
    review: Optional[ConditionReview] = None
    system: str = "condition_coworker.v1"

    @property
    def valid(self) -> bool:
        """Return the underlying chemistry/recommendation validity."""

        return self.result.valid

    def to_dict(self) -> Dict[str, Any]:
        """Serialize without truncating evidence, warnings, or provenance."""

        return {
            "system": self.system,
            "request": asdict(self.request),
            "result": self.result.to_dict(),
            "answer": self.answer,
            "review": self.review.to_dict() if self.review is not None else None,
        }


RetrosynthesisIssueCode = Literal[
    "functional_group_compatibility",
    "chemoselectivity",
    "ambiguous_reactive_site",
    "precursor_plausibility",
    "precursor_availability_unknown",
    "condition_feasibility",
    "protecting_group_requirement",
    "precedent_mismatch",
    "insufficient_evidence",
    "other",
]


@dataclass(frozen=True)
class RetrosynthesisRequest:
    """One chemistry-first, single-step retrosynthesis request."""

    target_smiles: str
    top_k: int = 5
    max_realizations_per_strategy: int = 3
    max_templates_to_apply: int = 500
    max_candidates_to_validate: int = 100
    use_context: bool = True
    include_l0: bool = True
    include_conditions: bool = True
    condition_top_k: int = 3
    condition_minimum_pool_size: Optional[int] = None
    unrestricted_condition_fallback: bool = False
    review: ConditionReviewSettings = field(default_factory=ConditionReviewSettings)

    def __post_init__(self) -> None:
        if not self.target_smiles.strip():
            raise ValueError("target_smiles must not be empty")
        if self.top_k < 1 or self.top_k > 50:
            raise ValueError("top_k must be between 1 and 50")
        if self.max_realizations_per_strategy < 1 or self.max_realizations_per_strategy > 10:
            raise ValueError("max_realizations_per_strategy must be between 1 and 10")
        if self.max_templates_to_apply < 1:
            raise ValueError("max_templates_to_apply must be positive")
        if self.max_candidates_to_validate < 1:
            raise ValueError("max_candidates_to_validate must be positive")
        if self.condition_top_k < 1 or self.condition_top_k > 10:
            raise ValueError("condition_top_k must be between 1 and 10")
        if (
            self.condition_minimum_pool_size is not None
            and self.condition_minimum_pool_size < 1
        ):
            raise ValueError("condition_minimum_pool_size must be positive")


@dataclass(frozen=True)
class RetrosynthesisStrategyCondition:
    """Condition evidence attached to one deterministic strategy."""

    strategy_id: str
    evidence: RetrosynthesisConditionEvidence

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the evidence without weakening its domain contract."""

        return {
            "strategy_id": self.strategy_id,
            "evidence": self.evidence.to_dict(),
        }


@dataclass(frozen=True)
class RetrosynthesisCandidateReview:
    """One evidence-linked advisory verdict for a validated strategy."""

    strategy_id: str
    original_rank: int
    suggested_rank: int
    verdict: ReviewVerdict
    issue_codes: Tuple[RetrosynthesisIssueCode, ...]
    evidence_ids: Tuple[str, ...]
    rationale: str
    confidence: float


@dataclass(frozen=True)
class RetrosynthesisReview:
    """Auditable outcome of bounded LLM review of validated strategies."""

    status: Literal["skipped", "completed", "failed"]
    provider: str
    model: str
    trigger_reasons: Tuple[str, ...] = ()
    summary: str = ""
    candidates: Tuple[RetrosynthesisCandidateReview, ...] = ()
    questions: Tuple[str, ...] = ()
    presentation_strategy_ids: Tuple[str, ...] = ()
    warning: Optional[str] = None
    response_id: Optional[str] = None
    input_tokens: int = 0
    output_tokens: int = 0
    provider_attempts: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serialize review evidence and provider metadata."""

        return asdict(self)


@dataclass(frozen=True)
class RetrosynthesisResponse:
    """Lossless deterministic strategies plus optional bounded review."""

    request: RetrosynthesisRequest
    valid: bool
    strategies: Tuple[StrategyProposal, ...] = ()
    condition_evidence: Tuple[RetrosynthesisStrategyCondition, ...] = ()
    answer: str = ""
    review: Optional[RetrosynthesisReview] = None
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    library_path: Optional[str] = None
    system: str = "chem_coworker.retrosynthesis.v1"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize all strategies, evidence, and model annotations."""

        return {
            "system": self.system,
            "request": asdict(self.request),
            "valid": self.valid,
            "strategies": [item.to_dict() for item in self.strategies],
            "condition_evidence": [
                item.to_dict() for item in self.condition_evidence
            ],
            "answer": self.answer,
            "review": self.review.to_dict() if self.review is not None else None,
            "warnings": list(self.warnings),
            "error": self.error,
            "library_path": self.library_path,
        }


MultistepIssueCode = Literal[
    "cross_step_functional_group_compatibility",
    "chemoselectivity",
    "protecting_group_strategy",
    "redox_or_protection_churn",
    "starting_material_evidence",
    "condition_feasibility",
    "route_convergence",
    "precedent_mismatch",
    "insufficient_evidence",
    "user_preference_mismatch",
    "other",
]


@dataclass(frozen=True)
class MultistepRetrosynthesisRequest:
    """One bounded, chemistry-first multistep route-search request."""

    target_smiles: str
    top_k: int = 5
    max_depth: Literal[2, 3] = 3
    molecular_weight_threshold: float = 150.0
    per_step_top_k: int = 5
    beam_width: int = 20
    max_expansions: int = 12
    max_templates_to_apply: int = 40
    max_candidates_to_validate: int = 10
    use_context: bool = True
    include_l0: bool = True
    include_conditions: bool = True
    condition_top_k: int = 3
    strategic_guidance: str = ""
    review: ConditionReviewSettings = field(default_factory=ConditionReviewSettings)

    def __post_init__(self) -> None:
        if not self.target_smiles.strip():
            raise ValueError("target_smiles must not be empty")
        if self.top_k < 1 or self.top_k > 10:
            raise ValueError("multistep top_k must be between 1 and 10")
        if self.max_depth not in {2, 3}:
            raise ValueError("multistep max_depth must be 2 or 3")
        if self.molecular_weight_threshold <= 0.0:
            raise ValueError("molecular_weight_threshold must be positive")
        for value, name in (
            (self.per_step_top_k, "per_step_top_k"),
            (self.beam_width, "beam_width"),
            (self.max_expansions, "max_expansions"),
            (self.max_templates_to_apply, "max_templates_to_apply"),
            (self.max_candidates_to_validate, "max_candidates_to_validate"),
            (self.condition_top_k, "condition_top_k"),
        ):
            if value < 1:
                raise ValueError(f"{name} must be positive")
        if self.condition_top_k > 10:
            raise ValueError("condition_top_k must be at most 10")
        if len(self.strategic_guidance) > 2_000:
            raise ValueError("strategic_guidance must be at most 2000 characters")


@dataclass(frozen=True)
class MultistepRouteReview:
    """One evidence-linked advisory verdict for an existing route."""

    route_id: str
    original_rank: int
    suggested_rank: int
    verdict: ReviewVerdict
    issue_codes: Tuple[MultistepIssueCode, ...]
    evidence_ids: Tuple[str, ...]
    rationale: str
    confidence: float


@dataclass(frozen=True)
class MultistepReview:
    """Auditable route-level LLM review that cannot mutate route chemistry."""

    status: Literal["skipped", "completed", "failed"]
    provider: str
    model: str
    reviewed_route_kind: Literal["solved", "partial", "none"] = "none"
    trigger_reasons: Tuple[str, ...] = ()
    summary: str = ""
    routes: Tuple[MultistepRouteReview, ...] = ()
    questions: Tuple[str, ...] = ()
    presentation_route_ids: Tuple[str, ...] = ()
    warning: Optional[str] = None
    response_id: Optional[str] = None
    input_tokens: int = 0
    output_tokens: int = 0
    provider_attempts: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Serialize review evidence and provider metadata."""

        return asdict(self)


@dataclass(frozen=True)
class MultistepRetrosynthesisResponse:
    """Deterministic route-search output plus optional route-level review."""

    request: MultistepRetrosynthesisRequest
    valid: bool
    result: Optional[MultistepRetrosynthesisResult] = None
    answer: str = ""
    review: Optional[MultistepReview] = None
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    library_path: Optional[str] = None
    stock_path: Optional[str] = None
    system: str = "chem_coworker.multistep_retrosynthesis.v1"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize route search and review without discarding evidence."""

        return {
            "system": self.system,
            "request": asdict(self.request),
            "valid": self.valid,
            "result": self.result.to_dict() if self.result is not None else None,
            "answer": self.answer,
            "review": self.review.to_dict() if self.review is not None else None,
            "warnings": list(self.warnings),
            "error": self.error,
            "library_path": self.library_path,
            "stock_path": self.stock_path,
        }
