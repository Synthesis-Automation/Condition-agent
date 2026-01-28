"""
Detection validation using reacted motifs patterns.

This module provides second-pass validation to correct common misclassifications
by using actual motif consumption patterns (reacted/formed motifs).
"""

from typing import Dict, Any, Set, List, Optional


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
    
    # Pattern-based validation rules for cross-coupling reactions
    
    # 1. SUZUKI-MIYAURA: Ar-B(OH)2 + Ar-Halide -> Ar-Ar
    if _has_organoboron(reacted_set) and _has_aryl_halide(reacted_set):
        if any(p in formed_set for p in ["Ar-Ar", "Ar-Alkenyl", "Alkenyl-Alkenyl"]):
            if initial_detection != "Suzuki_miyaura":
                return {
                    "reaction_type": "Suzuki_miyaura",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Organoboron + aryl halide → biaryl pattern",
                }
    
    # 2. BUCHWALD-HARTWIG (C-N Coupling): Ar-Halide + Amine -> Ar-NHR
    if _has_aryl_halide(reacted_set) and _has_amine(reacted_set):
        if any(p in formed_set for p in ["Ar-NHR", "Ar-NR2"]):
            if initial_detection not in ["Buchwald_CN", "C_N_Coupling"]:
                return {
                    "reaction_type": "C_N_Coupling",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Aryl halide + amine → arylamine pattern",
                }
    
    # 3. STILLE: Ar-Halide + Ar-SnR3 -> Ar-Ar
    if _has_organotin(reacted_set) and _has_aryl_halide(reacted_set):
        if "Ar-Ar" in formed_set or "Ar-Alkenyl" in formed_set:
            if initial_detection != "Stille":
                return {
                    "reaction_type": "Stille",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Organotin + aryl halide pattern",
                }
    
    # 4. NEGISHI: Ar-Halide + Ar-ZnX -> Ar-Ar
    if _has_organozinc(reacted_set) and _has_aryl_halide(reacted_set):
        if "Ar-Ar" in formed_set:
            if initial_detection != "Negishi":
                return {
                    "reaction_type": "Negishi",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Organozinc + aryl halide pattern",
                }
    
    # 5. HECK: Ar-Halide + Alkenyl -> Ar-Alkenyl
    if _has_aryl_halide(reacted_set) and "Alkenyl" in reacted_set:
        if "Ar-Alkenyl" in formed_set:
            if initial_detection != "Heck":
                return {
                    "reaction_type": "Heck",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Aryl halide + alkene pattern",
                }
    
    # 6. SONOGASHIRA: Ar-Halide + Alkynyl -> Ar-Alkynyl
    if _has_aryl_halide(reacted_set) and any(m in reacted_set for m in ["Alkynyl", "Ar-Alkynyl"]):
        if "Ar-Alkynyl" in formed_set:
            if initial_detection != "Sonogashira":
                return {
                    "reaction_type": "Sonogashira",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Aryl halide + alkyne pattern",
                }
    
    # Exclusion rule: Arylation_Ar_H with organometallic nucleophile
    if initial_detection == "Arylation_Ar_H":
        if _has_organometallic_nucleophile(reacted_set):
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


def _has_organoboron(motifs: Set[str]) -> bool:
    """Check for organoboron nucleophiles."""
    return any(m in motifs for m in ["Ar-B(OH)2", "Ar-B(OR)2", "Alkyl-B(OH)2", "Alkenyl-B(OH)2"])


def _has_aryl_halide(motifs: Set[str]) -> bool:
    """Check for aryl halide electrophiles."""
    return any(m in motifs for m in ["Ar-Br", "Ar-Cl", "Ar-I", "Ar-F"])


def _has_amine(motifs: Set[str]) -> bool:
    """Check for amine nucleophiles."""
    return any(m in motifs for m in ["Ar-NH2", "RNH2", "Ar-NHR", "R2NH"])


def _has_organotin(motifs: Set[str]) -> bool:
    """Check for organotin nucleophiles."""
    return any(m in motifs for m in ["Ar-SnR3", "Alkenyl-SnR3", "Alkyl-SnR3"])


def _has_organozinc(motifs: Set[str]) -> bool:
    """Check for organozinc nucleophiles."""
    return any(m in motifs for m in ["Ar-ZnX", "Alkyl-ZnX", "Alkenyl-ZnX"])


def _has_organometallic_nucleophile(motifs: Set[str]) -> bool:
    """Check for any organometallic nucleophile."""
    return (
        _has_organoboron(motifs)
        or _has_organotin(motifs)
        or _has_organozinc(motifs)
    )
