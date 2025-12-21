"""
Unified reaction detection API (taxonomy v2).

This module exposes ``detect_reaction`` for compatibility, but detection now
uses motif-slot evidence from ``chemtools.taxonomy.v2``. ML/SMARTS toggles are
accepted for backwards compatibility and ignored.
"""

from __future__ import annotations

from typing import Any, Dict, List

from .taxonomy.v2.reaction_detection import detect_reaction_types


def _confidence_from_match(matched_slots: int, required_slots: int) -> float:
    if required_slots <= 0:
        return 0.0
    return float(matched_slots) / float(required_slots)


def _build_result(
    *,
    matches: List[Dict[str, Any]],
    best: Dict[str, Any] | None,
    error: str | None,
) -> Dict[str, Any]:
    details: Dict[str, Any] = {
        "slot_evidence": best.get("slot_evidence", {}) if best else {},
        "matches": matches,
    }
    if best:
        details["matched_slots"] = best.get("matched_slots")
        details["required_slots"] = best.get("required_slots")
    if error:
        details["error"] = error

    family = best.get("reaction_type") if best else "Unknown"
    confidence = _confidence_from_match(
        int(best.get("matched_slots", 0)) if best else 0,
        int(best.get("required_slots", 0)) if best else 0,
    )

    return {
        "family": family or "Unknown",
        "confidence": confidence,
        "method": "motif_slots",
        "agreement": None,
        "status": "rule_only",
        "details": details,
    }


def detect_reaction(
    reaction_smiles: str,
    use_ml: bool = True,
    use_taxonomy_smarts: bool = True,
) -> Dict[str, Any]:
    """
    Detect reaction family from reaction SMILES using taxonomy v2.

    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products).
        use_ml: Ignored (kept for compatibility with old API).
        use_taxonomy_smarts: Ignored (kept for compatibility with old API).

    Returns:
        Dict with ``family``, ``confidence``, ``method``, and slot evidence
        metadata under ``details``.
    """
    _ = use_ml, use_taxonomy_smarts
    result = detect_reaction_types(reaction_smiles)
    matches = [match.to_dict() for match in result.matches]
    best = matches[0] if matches else None
    return _build_result(matches=matches, best=best, error=result.error)
