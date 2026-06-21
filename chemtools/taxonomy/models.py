"""Shared taxonomy model aliases and simple value objects."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

from .reaction_catalog import ReactionTypeDefinition, SlotRequirement


@dataclass(frozen=True)
class MotifDefinition:
    """Canonical motif entry exposed by taxonomy services."""

    motif_id: str
    payload: Mapping[str, Any]

    @property
    def label(self) -> str:
        return str(self.payload.get("label") or self.payload.get("name") or self.motif_id)

    @property
    def group_a(self) -> Optional[str]:
        value = self.payload.get("A") or self.payload.get("group_a")
        return str(value) if value else None

    @property
    def group_b(self) -> Optional[str]:
        value = self.payload.get("B") or self.payload.get("group_b")
        return str(value) if value else None


TaxonomyPayload = Dict[str, Any]

__all__ = [
    "MotifDefinition",
    "ReactionTypeDefinition",
    "SlotRequirement",
    "TaxonomyPayload",
]
