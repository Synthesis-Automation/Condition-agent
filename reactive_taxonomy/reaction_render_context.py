"""Shared terminal context for reaction text and graphic rendering."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from .notation import notation_style
from .reaction_models import (
    PartialProductTransformation,
    ReactionAnalysis,
    ReactionComponent,
    ReactionCoreProjection,
    ReactionInterpretation,
    ReactionObservation,
    ReactionSignature,
    ReactionTopology,
)


@dataclass(frozen=True)
class ReactionRenderContext:
    """Completed chemistry shared by sibling terminal renderers."""

    observation: ReactionObservation
    reactants: Tuple[ReactionComponent, ...]
    agents: Tuple[ReactionComponent, ...]
    products: Tuple[ReactionComponent, ...]
    signature: Optional[ReactionSignature]
    partial_transformation: Optional[PartialProductTransformation]
    interpretation: Optional[ReactionInterpretation]
    style: str
    schema_version: str = "1.0"

    @property
    def reaction_core(self) -> Optional[ReactionCoreProjection]:
        return self.observation.core

    @property
    def reaction_topology(self) -> Optional[ReactionTopology]:
        return self.observation.topology


def build_reaction_render_context(
    *,
    observation: ReactionObservation,
    reactants: Tuple[ReactionComponent, ...] = (),
    agents: Tuple[ReactionComponent, ...] = (),
    products: Tuple[ReactionComponent, ...] = (),
    signature: Optional[ReactionSignature] = None,
    partial_transformation: Optional[PartialProductTransformation] = None,
    interpretation: Optional[ReactionInterpretation] = None,
    style: str = "unicode",
) -> ReactionRenderContext:
    """Build the sole typed input contract for terminal reaction rendering."""
    notation_style(style)
    return ReactionRenderContext(
        observation=observation,
        reactants=reactants,
        agents=agents,
        products=products,
        signature=signature,
        partial_transformation=partial_transformation,
        interpretation=interpretation,
        style=style,
    )


def reaction_render_context_from_analysis(
    analysis: ReactionAnalysis,
    *,
    style: Optional[str] = None,
) -> ReactionRenderContext:
    """Project one completed analysis onto the shared rendering contract."""
    if analysis.observation is None:
        raise ValueError("Reaction rendering requires a ReactionObservation")
    selected_style = style or (
        analysis.reaction_label.style if analysis.reaction_label else "unicode"
    )
    return build_reaction_render_context(
        observation=analysis.observation,
        reactants=analysis.reactants,
        agents=analysis.agents,
        products=analysis.products,
        signature=analysis.reaction_signature,
        partial_transformation=analysis.partial_product_transformation,
        interpretation=analysis.interpretation,
        style=selected_style,
    )


__all__ = [
    "ReactionRenderContext",
    "build_reaction_render_context",
    "reaction_render_context_from_analysis",
]
