"""Public request and response contracts for the condition coworker."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Optional, Tuple

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
        }
