"""Public request and response contracts for the condition coworker."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Literal, Optional, Tuple

from condition_recommender import GenericRecommendationResult


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
    max_candidates: int = 5
    max_output_tokens: int = 2_000
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
class ConditionReview:
    """Auditable outcome of the optional LLM review stage."""

    status: Literal["skipped", "completed", "failed"]
    provider: str
    model: str
    trigger_reasons: Tuple[str, ...] = ()
    summary: str = ""
    candidates: Tuple[ConditionCandidateReview, ...] = ()
    questions: Tuple[str, ...] = ()
    presentation_recipe_ids: Tuple[str, ...] = ()
    warning: Optional[str] = None
    response_id: Optional[str] = None
    input_tokens: int = 0
    output_tokens: int = 0

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
