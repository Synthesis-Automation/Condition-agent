"""
Detection validation using reacted motifs patterns.

This module provides second-pass validation to correct common misclassifications
by using actual motif consumption patterns (reacted/formed motifs).

NOW TAXONOMY-DRIVEN: Validation patterns are loaded from taxonomy files
(reaction_types.v4.0.json + compound_logic.json) instead of hardcoded.
"""

from typing import Dict, Any, Set, List, Optional, Tuple
from chemtools.taxonomy.reaction_catalog import (
    load_reaction_catalog,
    ReactionTypeDefinition,
)


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
    reacted_set = set(reacted_motifs or [])
    formed_set = set(formed_motifs or [])
    
    # Load taxonomy
    definitions, _ = _get_catalog()
    
    # TAXONOMY-DRIVEN VALIDATION
    # Check each reaction definition in the catalog for pattern matches
    for reaction_id, defn in definitions.items():
        if not defn.reactants:
            continue
            
        # Check if all reactant slots match the reacted motifs
        all_slots_match = True
        matched_slots = []
        
        for slot_name, slot_req in defn.reactants.items():
            if not slot_req.allowed:
                continue  # Skip empty slots
                
            if _motifs_match_slot(reacted_set, slot_req.allowed):
                matched_slots.append(slot_name)
            else:
                all_slots_match = False
                break
        
        # If all required reactant slots matched, check product formation (if defined)
        if all_slots_match and len(matched_slots) >= 2:  # At least 2 reactant slots
            # Check if expected products are formed (if product patterns are defined)
            product_match = True
            if defn.products:
                product_match = False
                for slot_name, slot_req in defn.products.items():
                    if slot_req.allowed and _motifs_match_slot(formed_set, slot_req.allowed):
                        product_match = True
                        break
            
            # If pattern matches and it's different from initial detection, correct it
            if product_match and reaction_id != initial_detection:
                return {
                    "reaction_type": reaction_id,
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": f"Taxonomy pattern: {' + '.join(matched_slots)} → {defn.name}",
                }
    
    # Exclusion rule: Arylation_Ar_H with organometallic nucleophile
    if initial_detection == "Arylation_Ar_H":
        if _has_organometallic_nucleophile_taxonomy(reacted_set, definitions):
            return {
                "reaction_type": "Unknown",
                "confidence": 0.3,
                "validation_method": "reacted_motifs_exclusion",
                "corrected_from": initial_detection,
                "reason": "Organometallic nucleophile present - not C-H activation",
            }
    
    # No correction needed
    return {
        "reaction_type": initial_detection,
        "confidence": initial_confidence,
        "validation_method": "slot_based_confirmed",
        "corrected_from": None,
        "reason": "Pattern consistent with slot-based detection",
    }


def _has_organometallic_nucleophile_taxonomy(
    motifs: Set[str], 
    definitions: Dict[str, ReactionTypeDefinition]
) -> bool:
    """
    Check for any organometallic nucleophile using taxonomy definitions.
    Looks for motifs matching organoboron, organotin, organozinc patterns.
    """
    # Find organoboron, organotin, organozinc patterns from taxonomy
    organometallic_motifs = set()
    
    for reaction_id in ["Suzuki_miyaura", "Stille", "Negishi"]:
        if reaction_id in definitions:
            defn = definitions[reaction_id]
            for slot_name, slot_req in defn.reactants.items():
                if "nucleophile" in slot_name.lower() or "organo" in slot_name.lower():
                    organometallic_motifs.update(slot_req.allowed)
    
    return any(m in motifs for m in organometallic_motifs)
