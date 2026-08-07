"""Versioned HTTP request and response contracts for the local web app."""

from __future__ import annotations

from typing import Any, Dict, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator


API_SCHEMA_VERSION = "1.0"


class StrictRequest(BaseModel):
    """Reject misspelled or unsupported request fields."""

    model_config = ConfigDict(extra="forbid")


class PrepareReactionRequest(StrictRequest):
    """A reaction drawn or pasted by the user."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)


class CompletionChoiceRequest(StrictRequest):
    """One explicit user choice for a proposed missing fragment source."""

    requirement_id: str = Field(min_length=1, max_length=200)
    option_id: Optional[str] = Field(default=None, max_length=200)
    custom_identifier: Optional[str] = Field(default=None, max_length=500)

    @model_validator(mode="after")
    def exactly_one_source(self) -> "CompletionChoiceRequest":
        if bool(self.option_id) == bool(self.custom_identifier):
            raise ValueError(
                "provide exactly one of option_id or custom_identifier"
            )
        return self


class RankingPreferencesRequest(StrictRequest):
    """Chemist-selected ranking profile and optional custom priorities."""

    profile_id: str = Field(default="default", min_length=1, max_length=100)
    weights: Dict[str, float] = Field(default_factory=dict)

    @model_validator(mode="after")
    def positive_custom_weight(self) -> "RankingPreferencesRequest":
        if any(value < 0.0 for value in self.weights.values()):
            raise ValueError("ranking weights must be non-negative")
        if self.weights and sum(self.weights.values()) <= 0.0:
            raise ValueError("at least one ranking weight must be positive")
        return self


class RecommendationRequest(StrictRequest):
    """One condition-recommendation request."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)
    top_k: int = Field(default=5, ge=1, le=50)
    minimum_pool_size: Optional[int] = Field(default=None, ge=1, le=100)
    unrestricted_fallback: bool = False
    use_rxnmapper: bool = True
    ranking_preferences: RankingPreferencesRequest = Field(
        default_factory=RankingPreferencesRequest
    )
    completion_choices: tuple[CompletionChoiceRequest, ...] = ()


class DiscoveryRequest(StrictRequest):
    """One structurally related precedent-search request."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)
    top_k: int = Field(default=10, ge=1, le=50)
    view: Literal[
        "closest_chemistry",
        "diverse_strategies",
        "successful_precedents",
        "failure_informed",
    ] = "closest_chemistry"
    include_low_yield: bool = True
    include_unreported_outcomes: bool = True
    use_rxnmapper: bool = True
    include_review: bool = False


class RenderReactionRequest(StrictRequest):
    """Server-side RDKit rendering request."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)
    width: int = Field(default=760, ge=240, le=2400)
    height: int = Field(default=220, ge=100, le=1200)


class ApiEnvelope(BaseModel):
    """Stable wrapper around versioned domain result payloads."""

    api_schema_version: str = API_SCHEMA_VERSION
    data: Dict[str, Any]


def envelope(data: Dict[str, Any]) -> Dict[str, Any]:
    """Wrap a response without copying nested domain contracts."""

    return {
        "api_schema_version": API_SCHEMA_VERSION,
        "data": data,
    }


__all__ = [
    "API_SCHEMA_VERSION",
    "ApiEnvelope",
    "CompletionChoiceRequest",
    "DiscoveryRequest",
    "PrepareReactionRequest",
    "RankingPreferencesRequest",
    "RecommendationRequest",
    "RenderReactionRequest",
    "envelope",
]

