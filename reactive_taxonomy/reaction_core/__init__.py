"""Type-agnostic reaction-core observation API."""

from .builder import build_reaction_core_projection
from .models import (
    ReactionCoreAttachmentPort,
    ReactionCoreAtomState,
    ReactionCoreAtomTransition,
    ReactionCoreEvent,
    ReactionCoreProjection,
    ReactionCoreQuality,
    ReactionCoreRemoteSubgraph,
    ReactionCoreStateChange,
)

__all__ = [
    "ReactionCoreAttachmentPort",
    "ReactionCoreAtomState",
    "ReactionCoreAtomTransition",
    "ReactionCoreEvent",
    "ReactionCoreProjection",
    "ReactionCoreQuality",
    "ReactionCoreRemoteSubgraph",
    "ReactionCoreStateChange",
    "build_reaction_core_projection",
]
