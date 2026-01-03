"""
Context-aware reactant classification using reaction type information.

This module provides a two-pass approach for identifying reactants:
1. Pass 1: Detect reaction type from reaction SMILES (or accept user-provided type)
2. Pass 2: Classify reactants with knowledge of expected functional groups for that reaction

This solves the multi-functional group problem where a molecule like Brc1ccc(N)cc1
contains both ArBr and ArNH2 - we can determine which is reactive based on reaction context.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .reactants import ReactantMatch, iter_reactant_matches
from .smiles import normalize_reaction
from ._registry import get_registry
from ...reaction_type_detection import detect_reaction  # New unified API


@dataclass
class ContextualReactantMatch:
    """Reactant match with reaction context and role information."""
    
    # Core match info
    match: ReactantMatch
    
    # Context info
    role: Optional[str] = None  # "electrophile", "nucleophile", "coupling_partner", etc.
    is_expected: bool = False   # True if this matches the expected reactant type for this reaction
    position: int = 0           # Position in reactant list (0-indexed)
    
    # Alternative matches (other functional groups in same molecule)
    alternative_matches: List[ReactantMatch] = None
    
    def __post_init__(self):
        if self.alternative_matches is None:
            self.alternative_matches = []
    
    @property
    def category(self) -> str:
        """Convenience accessor for match category."""
        return self.match.category
    
    @property
    def member_type(self) -> str:
        """Convenience accessor for match member type."""
        return self.match.member_type
    
    @property
    def name(self) -> str:
        """Convenience accessor for match name."""
        return self.match.name


@dataclass
class ReactionClassification:
    """Complete reaction classification with context-aware reactant assignments."""
    
    # Reaction identification
    reaction_smiles: str
    reaction_type: str
    reaction_confidence: float
    
    # Reactant classifications
    reactants: List[ContextualReactantMatch]
    
    # Metadata
    detection_method: str  # "user_provided", "auto_detected", "rule_based"
    has_multi_functional: bool = False  # True if any reactant has multiple functional groups


def _infer_reactant_role(reactant_type_id: str, reaction_type: str) -> Optional[str]:
    """
    Infer likely role of a reactant type in a given reaction.
    
    This is a heuristic based on common reaction patterns.
    """
    # Common electrophile patterns
    if reactant_type_id in {
        "ArX*", "Alkyl-X", "VinylX*", "HetAr-X", "Allylic-X", "Benzylic-X",
        "Ar-X", "Ar-Br", "Ar-Cl", "Ar-I", "Ar-F",
        "R-X", "R-Br", "R-Cl", "R-I", "R-F",
        "Vinyl-X", "Vinyl-Br", "Vinyl-Cl", "Vinyl-I",
        "HetAr-X", "HetAr-Br", "HetAr-Cl", "HetAr-I",
        "Allyl-X", "Allyl-Br", "Allyl-Cl", "Allyl-I",
        "Bn-X", "Bn-Br", "Bn-Cl", "Bn-I",
        "Propargyl-X", "Propargyl-Br", "Propargyl-Cl", "Propargyl-I", "Propargyl-OTf", "Propargyl-Sulfonate",
        "Alkynyl-X", "Alkynyl-Br", "Alkynyl-Cl", "Alkynyl-I", "Alkynyl-OTf", "Alkynyl-Sulfonate",
        "Csp2-X", "Csp2-Sulfonate", "Any-Sulfonate", "Ar-OTf", "Ar-OTs", "Ar-OMs",
        "R-OTf", "R-OTs", "R-OMs"
    }:
        return "electrophile"
    
    # Common nucleophile patterns
    if reactant_type_id in {
        "RNH2/R2NH", "ArNH2/Ar2NH", "ROH", "ArOH", "RSH",
        "Any-NH2", "Any-NHR", "Ar-NH2", "Ar-NHR", "Ar-NHAr", "AromN-H",
        "R-OH", "Ar-OH", "R2CH-OH", "RCH2-OH", "Any-SH",
        "Bn-NH2", "Bn-NHR", "Allyl-NH2", "Allyl-NHR", "Propargyl-NH2", "Propargyl-NHR",
        "Bn-OH", "Allyl-OH", "Propargyl-OH"
    }:
        return "nucleophile"
    
    # Organometallic coupling partners
    if reactant_type_id in {
        "ArB*", "RB*", "RMgX", "RZnX", "RLi", "R-M",
        "Ar-B(OH)2", "Ar-B(OR)2", "R-B(OH)2", "R-B(OR)2", "Ar-BF3K", "R-BF3K",
        "Vinyl-B(OR)2", "R-MgX", "R-ZnX", "R-Li", "R-M", "Ar-MgX", "Ar-ZnX", "Ar-Li", "Ar-M"
    }:
        return "coupling_partner"
    
    # Carbonyl compounds (often electrophiles or undergo addition)
    if reactant_type_id in {
        "Aldehyde", "Ketone", "Any-Aldehyde", "Any-Ketone",
        "Any-CO2H", "Ar-CO2H", "R-CO2H", "Any-Ester", "Ar-Ester", "Any-CHO", "Any-CO-R"
    }:
        return "electrophile"
    
    # Alkenes (Heck, metathesis, addition)
    if reactant_type_id in {"Alkene", "Any-Alkene"}:
        if "heck" in reaction_type.lower():
            return "coupling_partner"
        return "substrate"
    
    # Alkynes (Sonogashira)
    if reactant_type_id in {"Alkyne", "Any-Alkyne"}:
        if "sonogashira" in reaction_type.lower():
            return "coupling_partner"
        return "substrate"
    
    return None


def classify_reactants_with_context(
    reaction_smiles: str,
    reaction_type: Optional[str] = None,
    auto_detect: bool = True
) -> ReactionClassification:
    """
    Classify reactants with reaction context (Two-Pass Approach).
    
    Args:
        reaction_smiles: Full reaction SMILES string
        reaction_type: Optional user-provided reaction type (e.g., "suzuki_miyaura")
                      If None and auto_detect=True, will detect automatically
        auto_detect: If True and reaction_type is None, auto-detect reaction type
    
    Returns:
        ReactionClassification with context-aware reactant assignments
    
    Example:
        >>> # User provides reaction type
        >>> result = classify_reactants_with_context(
        ...     "Brc1ccc(N)cc1.B(O)(O)c1ccccc1.[Pd]>>...",
        ...     reaction_type="suzuki_miyaura"
        ... )
        >>> result.reactants[0].category  # "ArX*" (ArBr is reactive in Suzuki)
        >>> result.reactants[0].role      # "electrophile"
        
        >>> # Auto-detect reaction type
        >>> result = classify_reactants_with_context(
        ...     "Brc1ccccc1.Nc1ccc(Br)cc1.[Pd]>>..."
        ... )
        >>> result.reaction_type          # "c_n_cross_coupling" (auto-detected)
        >>> result.reactants[1].category  # "ArNH2/Ar2NH" (NH2 is reactive in BH)
    """
    # Step 1: Normalize reaction
    normalized = normalize_reaction(reaction_smiles)
    reactant_smiles = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (normalized.get("reactants") or [])
    ]
    reactant_smiles = [s for s in reactant_smiles if s]
    
    # Step 2: Determine reaction type
    detection_method = "unknown"
    confidence = 0.0
    
    if reaction_type:
        # User provided
        detection_method = "user_provided"
        confidence = 1.0
    elif auto_detect:
        # Auto-detect using unified API
        detection_result = detect_reaction(reaction_smiles, use_ml=True)
        reaction_type = detection_result.get("family", "Unknown")
        confidence = detection_result.get("confidence", 0.0)
        
        # Determine detection method
        method = detection_result.get("method", "rule")
        if "ml" in method.lower():
            detection_method = "ml_detected"
        else:
            detection_method = "rule_based"
    else:
        # No reaction type available
        reaction_type = "Unknown"
        detection_method = "none"
        confidence = 0.0
    
    # Step 3: Get expected reactant types for this reaction
    registry = get_registry()
    expected_types: List[Dict[str, Any]] = []
    
    if registry and reaction_type != "Unknown":
        reaction_def = registry.reaction_types.get(reaction_type)
        if reaction_def:
            expected_types = [
                {
                    "reactant_type_id": r.reactant_type_id,
                    "notes": r.notes,
                    "tokens": r.original_tokens
                }
                for r in reaction_def.reactants
            ]
    
    # Step 4: Classify each reactant with context
    contextual_reactants: List[ContextualReactantMatch] = []
    has_multi_functional = False
    
    for i, smiles in enumerate(reactant_smiles):
        # Get ALL functional groups in this molecule
        all_matches = iter_reactant_matches(smiles)
        
        if not all_matches:
            continue
        
        # Check if multi-functional
        if len(all_matches) > 1:
            has_multi_functional = True
        
        # Try to find a match that's in the expected types
        best_match = None
        is_expected = False
        role = None
        alternatives = []
        
        if expected_types:
            # Match against expected reactant types
            for expected in expected_types:
                expected_id = expected["reactant_type_id"]
                
                for match in all_matches:
                    if match.category == expected_id:
                        best_match = match
                        is_expected = True
                        role = _infer_reactant_role(expected_id, reaction_type)
                        # Collect other matches as alternatives
                        alternatives = [m for m in all_matches if m != best_match]
                        break
                
                if best_match:
                    break
        
        # Fallback: use most specific match
        if not best_match and all_matches:
            # Sort by specificity (non-general first, then by SMARTS length)
            sorted_matches = sorted(
                all_matches,
                key=lambda m: (m.is_general, -m.specificity, m.member_type)
            )
            best_match = sorted_matches[0]
            role = _infer_reactant_role(best_match.category, reaction_type)
            alternatives = sorted_matches[1:]
        
        if best_match:
            contextual_reactants.append(
                ContextualReactantMatch(
                    match=best_match,
                    role=role,
                    is_expected=is_expected,
                    position=i,
                    alternative_matches=alternatives
                )
            )
    
    return ReactionClassification(
        reaction_smiles=reaction_smiles,
        reaction_type=reaction_type,
        reaction_confidence=confidence,
        reactants=contextual_reactants,
        detection_method=detection_method,
        has_multi_functional=has_multi_functional
    )


def get_reactant_summary(classification: ReactionClassification) -> Dict[str, Any]:
    """
    Generate a human-readable summary of the reaction classification.
    
    Args:
        classification: Result from classify_reactants_with_context
    
    Returns:
        Dict with summary information
    """
    summary = {
        "reaction_type": classification.reaction_type,
        "confidence": classification.reaction_confidence,
        "detection_method": classification.detection_method,
        "num_reactants": len(classification.reactants),
        "has_multi_functional_substrates": classification.has_multi_functional,
        "reactants": []
    }
    
    for r in classification.reactants:
        reactant_info = {
            "position": r.position,
            "category": r.category,
            "member_type": r.member_type,
            "name": r.name,
            "role": r.role,
            "is_expected": r.is_expected,
            "has_alternatives": len(r.alternative_matches) > 0
        }
        
        if r.alternative_matches:
            reactant_info["alternative_functional_groups"] = [
                {
                    "category": m.category,
                    "member_type": m.member_type,
                    "name": m.name
                }
                for m in r.alternative_matches
            ]
        
        summary["reactants"].append(reactant_info)
    
    return summary


__all__ = [
    "ContextualReactantMatch",
    "ReactionClassification",
    "classify_reactants_with_context",
    "get_reactant_summary",
]
