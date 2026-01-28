"""
Context-aware reactant classification using reaction type information.

This module provides a taxonomy-driven two-pass approach for identifying reactants:
1. Pass 1: Detect reaction type from reaction SMILES (or accept user-provided type)
2. Pass 2: Classify reactants with knowledge of expected functional groups for that reaction

This solves the multi-functional group problem where a molecule like Brc1ccc(N)cc1
contains both ArBr and ArNH2 - we can determine which is reactive based on reaction context.

All reactant roles (electrophile, nucleophile, coupling_partner, etc.) are defined in
the taxonomy (chemtools/taxonomy/data/reaction_types.v*.json), not hardcoded.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .reactants import ReactantMatch, iter_reactant_matches
from .smiles import normalize_reaction
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
    Infer likely role of a reactant type in a given reaction using taxonomy.
    
    This queries the reaction_types taxonomy to determine if the reactant_type_id
    matches any of the defined role slots (electrophile, nucleophile, etc.).
    
    Args:
        reactant_type_id: Reactant category ID (e.g., "Ar-Br", "Ar-B(OR)2")
        reaction_type: Reaction type ID (e.g., "suzuki_miyaura", "c_n_cross_coupling")
    
    Returns:
        Role string like "electrophile", "nucleophile", "coupling_partner", or None
    """
    from ...taxonomy import reaction_catalog
    
    # Load reaction catalog
    definitions, _ = reaction_catalog.load_reaction_catalog()
    
    # Normalize reaction type (try exact match, then case-insensitive)
    reaction_def = definitions.get(reaction_type)
    if not reaction_def:
        # Try case-insensitive lookup
        for rxn_id, rxn in definitions.items():
            if rxn_id.lower() == reaction_type.lower():
                reaction_def = rxn
                break
    
    if not reaction_def:
        return None
    
    # reaction_def.reactants is a dict: {"electrophile": SlotRequirement, "nucleophile": SlotRequirement, ...}
    # Each SlotRequirement has .allowed = list of reactant type IDs
    if not isinstance(reaction_def.reactants, dict):
        return None
    
    # Check each role slot to see if our reactant_type_id is in the allowed list
    for role, slot_req in reaction_def.reactants.items():
        if hasattr(slot_req, 'allowed'):
            if reactant_type_id in slot_req.allowed:
                return role
    
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
    
    # Step 3: Get expected reactant types for this reaction from taxonomy
    from ...taxonomy import reaction_catalog
    
    expected_types: List[str] = []  # List of allowed reactant type IDs
    
    if reaction_type != "Unknown":
        definitions, _ = reaction_catalog.load_reaction_catalog()
        
        # Try exact match, then case-insensitive
        reaction_def = definitions.get(reaction_type)
        if not reaction_def:
            for rxn_id, rxn in definitions.items():
                if rxn_id.lower() == reaction_type.lower():
                    reaction_def = rxn
                    break
        
        if reaction_def and isinstance(reaction_def.reactants, dict):
            # Collect all allowed reactant types from all roles
            for role, slot_req in reaction_def.reactants.items():
                if hasattr(slot_req, 'allowed'):
                    expected_types.extend(slot_req.allowed)
    
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
            # Match against expected reactant types from taxonomy
            for match in all_matches:
                if match.category in expected_types:
                    best_match = match
                    is_expected = True
                    role = _infer_reactant_role(match.category, reaction_type)
                    # Collect other matches as alternatives
                    alternatives = [m for m in all_matches if m != best_match]
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
