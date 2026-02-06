"""
Detection validation using reacted motifs patterns.

This module provides second-pass validation to correct common misclassifications
using taxonomy constraints from either:
1) reacted/formed motif sets, or
2) CRK reaction keys (preferred path).

NOW TAXONOMY-DRIVEN: Validation patterns are loaded from taxonomy files
(reaction_types.v4.0.json + compound_logic.json) instead of hardcoded.
"""

import re
from typing import Dict, Any, Set, List, Optional, Tuple
from chemtools.taxonomy.reaction_catalog import (
    load_reaction_catalog,
    ReactionTypeDefinition,
)
from .utils import normalize_motif_id


# Cache the taxonomy for performance
_CATALOG_CACHE: Optional[Tuple[Dict[str, ReactionTypeDefinition], Dict[str, str]]] = None


def _get_catalog() -> Tuple[Dict[str, ReactionTypeDefinition], Dict[str, str]]:
    """Load and cache the reaction catalog."""
    global _CATALOG_CACHE
    if _CATALOG_CACHE is None:
        _CATALOG_CACHE = load_reaction_catalog()
    return _CATALOG_CACHE


def _motifs_match_slot(motifs: Set[str], allowed: List[str]) -> bool:
    """Check if any motif matches the allowed patterns for a slot."""
    return any(m in motifs for m in allowed)


def _split_motif_tokens(value: str) -> List[str]:
    tokens = [tok.strip() for tok in re.split(r"[|,;/]", value) if tok.strip()]
    return [normalize_motif_id(tok) for tok in tokens if tok and tok != "[]"]


def _parse_crk_key(
    reaction_key: str,
) -> Tuple[Set[str], Set[str], Set[str], List[str], List[str]]:
    """
    Parse CRK key into reacted/formed/spectator motif sets and bond tokens.

    Returns:
        (reacted_set, formed_set, spectators_set, bond_formed_tokens, bond_broken_tokens)
    """
    if not reaction_key:
        return set(), set(), set(), [], []

    text = str(reaction_key).strip()
    if not text or "->" not in text:
        return set(), set(), set(), [], []

    sections = [section.strip() for section in text.split(" | ") if section.strip()]
    if not sections:
        return set(), set(), set(), [], []

    summary = sections[0]
    if summary.startswith("CRK-v1"):
        summary = summary[len("CRK-v1"):].strip()
    if summary.startswith("|"):
        summary = summary[1:].strip()

    reacted: Set[str] = set()
    formed: Set[str] = set()
    if "->" in summary:
        reactant_part, product_part = summary.split("->", 1)
        reacted = set(_split_motif_tokens(reactant_part.strip()))
        formed = set(_split_motif_tokens(product_part.strip()))

    spectators: Set[str] = set()
    formed_bonds: List[str] = []
    broken_bonds: List[str] = []
    for section in sections[1:]:
        lower = section.lower()
        if lower.startswith("spectators:"):
            payload = section.split(":", 1)[1].strip()
            spectators = set(_split_motif_tokens(payload))
        elif lower.startswith("bond_formed:"):
            payload = section.split(":", 1)[1].strip()
            formed_bonds = [tok.strip() for tok in payload.split(";") if tok.strip()]
        elif lower.startswith("bond_broken:"):
            payload = section.split(":", 1)[1].strip()
            broken_bonds = [tok.strip() for tok in payload.split(";") if tok.strip()]

    return reacted, formed, spectators, formed_bonds, broken_bonds


def _match_reaction_catalog(
    reacted_set: Set[str],
    formed_set: Set[str],
) -> Optional[Tuple[str, List[str], str]]:
    """Return best taxonomy match as (reaction_id, matched_slots, display_name)."""
    definitions, _ = _get_catalog()
    for reaction_id, defn in definitions.items():
        if not defn.reactants:
            continue

        all_slots_match = True
        matched_slots: List[str] = []
        for slot_name, slot_req in defn.reactants.items():
            if not slot_req.allowed:
                continue
            if _motifs_match_slot(reacted_set, slot_req.allowed):
                matched_slots.append(slot_name)
            else:
                all_slots_match = False
                break

        if not all_slots_match or len(matched_slots) < 1:
            continue

        product_match = True
        if defn.products:
            product_match = False
            for slot_req in defn.products.values():
                if slot_req.allowed and _motifs_match_slot(formed_set, slot_req.allowed):
                    product_match = True
                    break

        if product_match:
            return reaction_id, matched_slots, defn.name
    return None


def _format_validation_response(
    *,
    initial_detection: str,
    initial_confidence: float,
    match: Optional[Tuple[str, List[str], str]],
    method: str,
    reason_prefix: str,
) -> Dict[str, Any]:
    if match is None:
        return {
            "reaction_type": initial_detection,
            "confidence": initial_confidence,
            "validation_method": "slot_based_confirmed",
            "corrected_from": None,
            "reason": "Pattern consistent with slot-based detection",
        }

    reaction_id, matched_slots, display_name = match
    if reaction_id != initial_detection:
        return {
            "reaction_type": reaction_id,
            "confidence": 0.95,
            "validation_method": method,
            "corrected_from": initial_detection,
            "reason": f"{reason_prefix}: {' + '.join(matched_slots)} -> {display_name}",
        }

    return {
        "reaction_type": initial_detection,
        "confidence": max(initial_confidence, 0.95),
        "validation_method": "slot_based_confirmed",
        "corrected_from": None,
        "reason": "Pattern consistent with slot-based detection",
    }


def validate_detection_with_crk_key(
    initial_detection: str,
    initial_confidence: float,
    reaction_key: str,
) -> Dict[str, Any]:
    """
    Validate and refine reaction type detection from CRK key.

    This is the preferred streamlined path: CRK_raw -> taxonomy match.
    """
    reacted_set, formed_set, _spectators, _formed_bonds, _broken_bonds = _parse_crk_key(reaction_key)
    match = _match_reaction_catalog(reacted_set, formed_set)
    return _format_validation_response(
        initial_detection=initial_detection,
        initial_confidence=initial_confidence,
        match=match,
        method="crk_pattern",
        reason_prefix="Taxonomy pattern (CRK)",
    )


def validate_detection_with_reacted_motifs(
    initial_detection: str,
    initial_confidence: float,
    reacted_motifs: List[str],
    formed_motifs: List[str],
    spectator_motifs: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Validate and refine reaction type detection using reacted motifs.
    
    This adds a logical second-pass that corrects common misclassifications
    by using actual motif consumption patterns.
    
    NOW TAXONOMY-DRIVEN: Loads patterns from reaction_types.v4.0.json + compound_logic.json
    instead of hardcoding motif lists.
    
    Args:
        initial_detection: Reaction type from slot-based detection
        initial_confidence: Confidence from slot-based detection
        reacted_motifs: Motifs consumed in reaction (from aggregates)
        formed_motifs: Motifs formed in products (from aggregates)
        spectator_motifs: Motifs present but unchanged (optional)
        
    Returns:
        Dict with validated reaction_type, confidence, and validation metadata
    """
    reacted_set = {normalize_motif_id(str(m)) for m in (reacted_motifs or []) if m}
    formed_set = {normalize_motif_id(str(m)) for m in (formed_motifs or []) if m}
    match = _match_reaction_catalog(reacted_set, formed_set)
    return _format_validation_response(
        initial_detection=initial_detection,
        initial_confidence=initial_confidence,
        match=match,
        method="reacted_motifs_pattern",
        reason_prefix="Taxonomy pattern",
    )
