"""Thin tool wrappers exposed via the Condition MCP façade.

Each wrapper is intentionally deterministic and returns JSON-serialisable
payloads with an explicit ``schema_version`` field so they can be surfaced
through Model Context Protocol tool manifests.
"""

from .normalize import normalize_reaction
from .detect import detect_family
from .featurize import featurize_substrates

__all__ = [
    "normalize_reaction",
    "detect_family",
    "featurize_substrates",
]
