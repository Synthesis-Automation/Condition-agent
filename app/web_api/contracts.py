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
    library_mode: Literal["full", "compact"] = "full"
    top_k: int = Field(default=5, ge=1, le=50)
    minimum_pool_size: Optional[int] = Field(default=None, ge=1, le=100)
    unrestricted_fallback: bool = False
    use_rxnmapper: bool = True
    ranking_preferences: RankingPreferencesRequest = Field(
        default_factory=RankingPreferencesRequest
    )
    completion_choices: tuple[CompletionChoiceRequest, ...] = ()


class FeatureAnalysisRequest(StrictRequest):
    """One auto-detected molecule or reaction featurization request."""

    input_smiles: str = Field(min_length=1, max_length=20_000)
    use_rxnmapper: bool = True
    force_resolved_mapping: bool = False

    @model_validator(mode="after")
    def force_requires_mapping(self) -> "FeatureAnalysisRequest":
        if self.force_resolved_mapping and not self.use_rxnmapper:
            raise ValueError("forced resolved mapping requires RXNMapper")
        return self


class RetrosynthesisRequest(StrictRequest):
    """One structure-derived, single-step retrosynthesis request."""

    target_smiles: str = Field(min_length=1, max_length=20_000)
    library_mode: Literal["full", "compact"] = "compact"
    top_k: int = Field(default=10, ge=1, le=50)
    include_l0: bool = True
    use_context: bool = True
    diversify: bool = True
    use_precursor_realism: bool = False
    use_forward_validation: bool = True


class ForwardConditionProfileRequest(StrictRequest):
    """Coarse condition intent without pretending to be a resolved recipe."""

    strategy: Literal[
        "unspecified",
        "transition_metal_catalysis",
        "metal_free_polar",
        "radical",
        "photochemical",
        "thermal",
    ] = "unspecified"
    redox_mode: Literal["neutral", "oxidative", "reductive"] = "neutral"
    medium: Literal["neutral", "acidic", "basic"] = "neutral"
    catalyst_family: Literal[
        "unspecified",
        "palladium",
        "nickel",
        "copper",
        "iron",
        "cobalt",
        "rhodium",
        "ruthenium",
        "iridium",
        "other_transition_metal",
    ] = "unspecified"

    @model_validator(mode="after")
    def catalyst_requires_metal_strategy(self) -> "ForwardConditionProfileRequest":
        if (
            self.catalyst_family != "unspecified"
            and self.strategy != "transition_metal_catalysis"
        ):
            raise ValueError(
                "catalyst_family requires transition_metal_catalysis strategy"
            )
        return self


class ForwardSynthesisRequest(StrictRequest):
    """Predict products or audit one proposed forward route step."""

    starting_materials: str = Field(min_length=1, max_length=20_000)
    intended_product: Optional[str] = Field(default=None, max_length=20_000)
    operator_hint: Optional[str] = Field(default=None, max_length=500)
    recipe: Optional[Dict[str, Any]] = None
    condition_profile: Optional[ForwardConditionProfileRequest] = None
    library_mode: Literal["full", "compact"] = "compact"
    top_k: int = Field(default=10, ge=1, le=50)
    include_l0: bool = True
    include_self_reactions: bool = True

    @model_validator(mode="after")
    def molecule_collections_only(self) -> "ForwardSynthesisRequest":
        if ">" in self.starting_materials:
            raise ValueError(
                "starting_materials must be molecules, not reaction SMILES"
            )
        if self.intended_product and ">" in self.intended_product:
            raise ValueError(
                "intended_product must be one molecule, not reaction SMILES"
            )
        return self


class MultistepRetrosynthesisRequest(StrictRequest):
    """One bounded, structure-derived multistep retrosynthesis request."""

    target_smiles: str = Field(min_length=1, max_length=20_000)
    library_mode: Literal["full", "compact"] = "compact"
    top_k_routes: int = Field(default=5, ge=1, le=10)
    max_depth: Literal[2, 3] = 3
    molecular_weight_threshold: float = Field(default=150.0, gt=0, le=500)
    include_l0: bool = True
    use_context: bool = True
    diversify: bool = True
    use_precursor_realism: bool = False
    use_condition_availability: bool = False
    use_forward_validation: bool = True


class RetrosynthesisConditionsRequest(StrictRequest):
    """Progressive condition lookup for one retrosynthesis hit."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)
    library_mode: Literal["full", "compact"] = "compact"
    top_k: int = Field(default=3, ge=1, le=5)
    starting_materials: Optional[str] = Field(default=None, max_length=20_000)
    intended_product: Optional[str] = Field(default=None, max_length=20_000)
    operator_hint: Optional[str] = Field(default=None, max_length=500)
    use_forward_validation: bool = False
    include_l0: bool = True

    @model_validator(mode="after")
    def forward_audit_requires_step(self) -> "RetrosynthesisConditionsRequest":
        if self.use_forward_validation and not (
            self.starting_materials and self.intended_product
        ):
            raise ValueError(
                "forward validation requires starting_materials and intended_product"
            )
        return self


class RenderReactionRequest(StrictRequest):
    """Server-side RDKit rendering request."""

    reaction_smiles: str = Field(min_length=1, max_length=20_000)
    width: int = Field(default=760, ge=240, le=2400)
    height: int = Field(default=220, ge=100, le=1200)


class RenderMoleculeRequest(StrictRequest):
    """Server-side molecule rendering request."""

    molecule_smiles: str = Field(min_length=1, max_length=20_000)
    width: int = Field(default=760, ge=160, le=2400)
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
    "FeatureAnalysisRequest",
    "ForwardConditionProfileRequest",
    "ForwardSynthesisRequest",
    "MultistepRetrosynthesisRequest",
    "PrepareReactionRequest",
    "RankingPreferencesRequest",
    "RecommendationRequest",
    "RetrosynthesisRequest",
    "RetrosynthesisConditionsRequest",
    "RenderMoleculeRequest",
    "RenderReactionRequest",
    "envelope",
]
