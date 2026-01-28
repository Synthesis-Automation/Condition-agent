"""
Data structures and models for reaction type detection results.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class ReactionMatch:
    """A matched reaction type with evidence and confidence scoring."""
    
    reaction_type: str
    name: str
    category: Optional[str]
    slot_evidence: Dict[str, List[str]]
    matched_slots: int
    required_slots: int

    @property
    def confidence(self) -> float:
        """Confidence score as ratio of matched to required slots."""
        if self.required_slots == 0:
            return 0.0
        return self.matched_slots / self.required_slots

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        return {
            "reaction_type": self.reaction_type,
            "name": self.name,
            "category": self.category,
            "confidence": self.confidence,
            "matched_slots": self.matched_slots,
            "required_slots": self.required_slots,
            "slot_evidence": {slot: list(values) for slot, values in self.slot_evidence.items()},
        }


@dataclass(frozen=True)
class ReactionDetectionResult:
    """Overall detection result containing all matches and optional error."""
    
    matches: List[ReactionMatch]
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for serialization."""
        payload = {"matches": [match.to_dict() for match in self.matches]}
        if self.error:
            payload["error"] = self.error
        return payload
