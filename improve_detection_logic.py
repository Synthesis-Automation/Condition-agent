"""
Improved Reaction Detection Using Two-Pass Validation

This implements a logical improvement to reaction type detection by adding
a second validation pass that uses reacted_motifs information.
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from typing import Dict, Any, Set, List, Optional


def validate_detection_with_reacted_motifs(
    initial_detection: str,
    initial_confidence: float,
    reacted_motifs: List[str],
    formed_motifs: List[str],
    spectator_motifs: List[str],
) -> Dict[str, Any]:
    """
    Validate and refine reaction type detection using reacted motifs.
    
    This is a logical second-pass that corrects common misclassifications
    by using actual motif consumption patterns.
    
    Args:
        initial_detection: Reaction type from slot-based detection
        initial_confidence: Confidence from slot-based detection
        reacted_motifs: Motifs consumed in reaction (from aggregates)
        formed_motifs: Motifs formed in products (from aggregates)
        spectator_motifs: Motifs present but unchanged
        
    Returns:
        Dict with validated reaction_type, confidence, and validation_notes
    """
    reacted_set = set(reacted_motifs)
    formed_set = set(formed_motifs)
    
    # Pattern-based validation rules
    # These are high-confidence patterns that override slot-based detection
    
    # 1. SUZUKI-MIYAURA: Ar-B(OH)2 + Ar-Halide -> Ar-Ar
    if _has_pattern(reacted_set, ["Ar-B(OH)2", "Ar-B(OR)2"], ["Ar-Br", "Ar-Cl", "Ar-I"]):
        if "Ar-Ar" in formed_set or "Ar-Alkenyl" in formed_set:
            if initial_detection != "Suzuki_miyaura":
                return {
                    "reaction_type": "Suzuki_miyaura",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Organoboron + aryl electrophile pattern detected",
                }
    
    # 2. BUCHWALD-HARTWIG: Ar-Halide + Ar-NH2 -> Ar-NHR
    if _has_pattern(reacted_set, ["Ar-Br", "Ar-Cl", "Ar-I"], ["Ar-NH2", "Ar-NHR", "RNH2"]):
        if any(m in formed_set for m in ["Ar-NHR", "Ar-NR2"]):
            if initial_detection not in ["Buchwald_CN", "C_N_Coupling"]:
                return {
                    "reaction_type": "C_N_Coupling",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Aryl halide + amine coupling pattern detected",
                }
    
    # 3. STILLE: Ar-Halide + Ar-SnR3 -> Ar-Ar
    if _has_pattern(reacted_set, ["Ar-Br", "Ar-Cl", "Ar-I"], ["Ar-SnR3", "Alkenyl-SnR3"]):
        if "Ar-Ar" in formed_set:
            if initial_detection != "Stille":
                return {
                    "reaction_type": "Stille",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Organotin + aryl electrophile pattern detected",
                }
    
    # 4. NEGISHI: Ar-Halide + Ar-ZnX -> Ar-Ar
    if _has_pattern(reacted_set, ["Ar-Br", "Ar-Cl", "Ar-I"], ["Ar-ZnX", "Alkyl-ZnX"]):
        if "Ar-Ar" in formed_set:
            if initial_detection != "Negishi":
                return {
                    "reaction_type": "Negishi",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Organozinc + aryl electrophile pattern detected",
                }
    
    # 5. HECK: Ar-Halide + Alkenyl -> Ar-Alkenyl
    if _has_pattern(reacted_set, ["Ar-Br", "Ar-Cl", "Ar-I"], ["Alkenyl"]):
        if "Ar-Alkenyl" in formed_set:
            if initial_detection != "Heck":
                return {
                    "reaction_type": "Heck",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Aryl halide + alkene pattern detected",
                }
    
    # 6. SONOGASHIRA: Ar-Halide + Alkynyl -> Ar-Alkynyl
    if _has_pattern(reacted_set, ["Ar-Br", "Ar-Cl", "Ar-I"], ["Alkynyl", "Ar-Alkynyl"]):
        if "Ar-Alkynyl" in formed_set:
            if initial_detection != "Sonogashira":
                return {
                    "reaction_type": "Sonogashira",
                    "confidence": 0.95,
                    "validation_method": "reacted_motifs_pattern",
                    "corrected_from": initial_detection,
                    "reason": "Aryl halide + alkyne pattern detected",
                }
    
    # 7. SNAr: Ar-Halide (activated) + Nu -> Ar-Nu
    # Only validate if Ar-H is detected but shouldn't be (common false positive)
    if initial_detection == "Arylation_Ar_H":
        # Check if there's an organo-nucleophile that's NOT just Ar-H
        if any(m in reacted_set for m in ["Ar-B(OH)2", "Ar-B(OR)2", "Ar-SnR3", "Ar-ZnX", "RNH2", "Ar-NH2"]):
            # This is a cross-coupling, not C-H activation
            return {
                "reaction_type": "Unknown",  # Let role classification handle it
                "confidence": 0.5,
                "validation_method": "reacted_motifs_exclusion",
                "corrected_from": initial_detection,
                "reason": "Detected organo-nucleophile inconsistent with C-H activation",
            }
    
    # No correction needed - initial detection is consistent
    return {
        "reaction_type": initial_detection,
        "confidence": initial_confidence,
        "validation_method": "slot_based_confirmed",
        "corrected_from": None,
        "reason": "Reacted motifs pattern consistent with slot-based detection",
    }


def _has_pattern(
    reacted_set: Set[str],
    group1_options: List[str],
    group2_options: List[str],
) -> bool:
    """
    Check if reacted_set contains at least one motif from each group.
    
    This implements logical AND across groups, OR within groups:
    (group1[0] OR group1[1] OR ...) AND (group2[0] OR group2[1] OR ...)
    """
    has_group1 = any(opt in reacted_set for opt in group1_options)
    has_group2 = any(opt in reacted_set for opt in group2_options)
    return has_group1 and has_group2


# Integration example
if __name__ == "__main__":
    print("=" * 80)
    print("TESTING IMPROVED DETECTION LOGIC")
    print("=" * 80)
    
    # Test Case 1: Suzuki misclassified as Arylation
    result = validate_detection_with_reacted_motifs(
        initial_detection="Arylation_Ar_H",
        initial_confidence=1.0,
        reacted_motifs=["Ar-B(OH)2", "Ar-Br"],
        formed_motifs=["Ar-Ar"],
        spectator_motifs=["Pyridine"],
    )
    
    print("\nTest 1: Suzuki Misclassification")
    print(f"Initial: {result.get('corrected_from', 'N/A')}")
    print(f"Validated: {result['reaction_type']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Reason: {result['reason']}")
    
    # Test Case 2: Buchwald-Hartwig
    result = validate_detection_with_reacted_motifs(
        initial_detection="Arylation_Ar_H",
        initial_confidence=1.0,
        reacted_motifs=["Ar-Br", "Ar-NH2"],
        formed_motifs=["Ar-NHR"],
        spectator_motifs=[],
    )
    
    print("\nTest 2: Buchwald-Hartwig Correction")
    print(f"Initial: {result.get('corrected_from', 'N/A')}")
    print(f"Validated: {result['reaction_type']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Reason: {result['reason']}")
    
    # Test Case 3: Legitimate Arylation (no correction needed)
    result = validate_detection_with_reacted_motifs(
        initial_detection="Arylation_Ar_H",
        initial_confidence=0.9,
        reacted_motifs=["Ar-Br", "Ar-H"],
        formed_motifs=["Ar-Ar"],
        spectator_motifs=[],
    )
    
    print("\nTest 3: Legitimate Arylation (No Correction)")
    print(f"Initial: {result.get('corrected_from', 'N/A')}")
    print(f"Validated: {result['reaction_type']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Reason: {result['reason']}")
    
    print("\n" + "=" * 80)
    print("INTEGRATION APPROACH")
    print("=" * 80)
    print("""
This function should be called in featurize_reaction() after aggregates
are calculated:

    # Current flow
    detection = detect_reaction_types(...)
    reaction_type = format_reaction_type_summary(detection)
    
    # ... featurize reactants/products ...
    
    aggregates = aggregate_reaction_features(...)
    
    # NEW: Validate with reacted motifs
    validated = validate_detection_with_reacted_motifs(
        initial_detection=reaction_type,
        initial_confidence=detection.confidence,
        reacted_motifs=aggregates.get('reacted_motifs', []),
        formed_motifs=aggregates.get('formed_motifs', []),
        spectator_motifs=aggregates.get('spectator_motifs', []),
    )
    
    # Use validated result
    reaction_type = validated['reaction_type']
    confidence = validated['confidence']
    
    # Store validation info in metadata
    if validated.get('corrected_from'):
        reaction['validation'] = {
            'original_detection': validated['corrected_from'],
            'validation_method': validated['validation_method'],
            'validation_reason': validated['reason'],
        }
    
This approach:
✅ Doesn't break existing slot-based detection
✅ Adds logical validation layer using reaction patterns
✅ Provides clear audit trail (corrected_from, reason)
✅ High confidence corrections (0.95) for known patterns
✅ Extensible - easy to add more patterns
""")
