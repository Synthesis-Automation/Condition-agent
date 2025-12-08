"""Reaction family detection wrapper for Condition MCP."""

from __future__ import annotations

from typing import Any, Dict, List

from pydantic import BaseModel, Field

from chemtools.reaction_type_detection import detect_reaction  # New unified API

from .base import SchemaStamped, pick_first, validate_payload


def _normalise_reactants(reactants: List[str]) -> List[str]:
    out: List[str] = []
    for item in reactants:
        if isinstance(item, str) and item.strip():
            out.append(item.strip())
    return out


class DetectFamilyInput(BaseModel):
    """Input payload for the ``detect_family`` tool."""

    reactants: List[str] = Field(..., description="Canonical reactant SMILES strings")


class DetectFamilyOutput(SchemaStamped):
    """Detected reaction family with confidence metrics."""

    family: str
    confidence: float
    hits: Dict[str, bool]
    has_conflict: bool


def detect_family(data: Dict[str, Any]) -> Dict[str, Any]:
    """Detect a reaction family for a set of reactant SMILES strings."""

    payload = validate_payload(DetectFamilyInput, data)
    reactants = _normalise_reactants(payload.reactants)
    
    # Convert reactants to pseudo-reaction for unified API
    pseudo_reaction = ".".join(reactants) + ">>"
    result = detect_reaction(pseudo_reaction, use_ml=False)

    family = pick_first([result.get("family"), "Unknown"]) or "Unknown"
    confidence_raw = result.get("confidence", 0.0)
    try:
        confidence = float(confidence_raw)
    except Exception:  # pragma: no cover - defensive conversion
        confidence = 0.0
    
    # Extract hits from details
    hits = result.get("details", {}).get("functional_groups", {})

    output = DetectFamilyOutput(
        family=family,
        confidence=confidence,
        hits={str(k): bool(v) for k, v in hits.items()},
        has_conflict=bool(result.get("status") == "conflict"),
    )
    return output.model_dump()
