"""Reaction domain model aliases."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional

from .featurize import CrkResult
from .inference import GeneralReactionAnalysis, ReactionDecision, ReactionValidation
from .typing import DetectionResult, ReactionMatch


@dataclass(frozen=True)
class ReactionEvents:
    """Structured wrapper for reaction event evidence."""

    payload: Mapping[str, Any] = field(default_factory=dict)

    def to_payload(self) -> Dict[str, Any]:
        return dict(self.payload)


@dataclass(frozen=True)
class ReactionFeatures:
    """Structured wrapper for the current reaction feature bundle."""

    reaction_smiles: str
    payload: Mapping[str, Any] = field(default_factory=dict)
    reaction_type: Optional[str] = None
    reaction_key: Optional[str] = None

    @classmethod
    def from_payload(cls, payload: Mapping[str, Any]) -> "ReactionFeatures":
        reaction_type = payload.get("reaction_type")
        if isinstance(reaction_type, Mapping):
            reaction_type_value = reaction_type.get("reaction_type")
        else:
            reaction_type_value = reaction_type
        return cls(
            reaction_smiles=str(payload.get("reaction_smiles") or ""),
            payload=payload,
            reaction_type=str(reaction_type_value) if reaction_type_value else None,
            reaction_key=str(payload.get("reaction_key") or "") or None,
        )

    def to_payload(self) -> Dict[str, Any]:
        return dict(self.payload)

__all__ = [
    "CrkResult",
    "DetectionResult",
    "GeneralReactionAnalysis",
    "ReactionEvents",
    "ReactionFeatures",
    "ReactionDecision",
    "ReactionMatch",
    "ReactionValidation",
]
