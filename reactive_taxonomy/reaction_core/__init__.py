"""Type-agnostic reaction-core observation and presentation API.

The package keeps graph-derived identity separate from chemist-facing
rendering.  Callers should construct cores through
``build_reaction_core_projection`` and treat labels as presentation only.
"""

from .builder import build_reaction_core_projection
from .models import (
    ReactionCoreAbstraction,
    ReactionCoreAttachmentPort,
    ReactionCoreAtomState,
    ReactionCoreAtomTransition,
    ReactionCoreEvent,
    ReactionCoreProjection,
    ReactionCoreRemoteSubgraph,
)

__all__ = [
    "ReactionCoreAbstraction",
    "ReactionCoreAttachmentPort",
    "ReactionCoreAtomState",
    "ReactionCoreAtomTransition",
    "ReactionCoreEvent",
    "ReactionCoreProjection",
    "ReactionCoreRemoteSubgraph",
    "build_reaction_core_projection",
]
