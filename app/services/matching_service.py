"""
Matching Service - SMILES normalization and reaction family detection.

This service handles:
- SMILES normalization (canonicalization)
- Reaction family detection using router
- Reaction type detection using rxn-insight + router fallback
"""

from typing import Dict, Any, List
from chemtools import chem
from chemtools.contracts import NormalizeRequest, DetectFamilyRequest, DetectTypeRequest
from chemtools.exceptions import ValidationError

try:
    from chemtools.reaction_type_detector import (
        detect_reaction_type as rxn_detect_type,
        is_available as rxn_insight_available,
    )
    _HAS_RXN_INSIGHT = True
except Exception:
    _HAS_RXN_INSIGHT = False


def normalize_smiles(req: NormalizeRequest) -> str:
    """
    Normalize a SMILES string to canonical form.
    
    Args:
        req: NormalizeRequest with SMILES string
        
    Returns:
        Canonical SMILES string
        
    Raises:
        ValidationError: If SMILES is invalid
    """
    if not req.smiles or not req.smiles.strip():
        raise ValidationError("SMILES cannot be empty")
    
    return chem.smiles.normalize(req.smiles)


def detect_family(req: DetectFamilyRequest) -> Dict[str, Any]:
    """
    Detect reaction family from reactants.
    
    Args:
        req: DetectFamilyRequest with reactants list
        
    Returns:
        Dictionary with family detection results
        
    Raises:
        ValidationError: If reactants are invalid
    """
    if not req.reactants:
        raise ValidationError("Reactants list cannot be empty")
    
    return chem.router.detect_family(req.reactants)


def detect_reaction_type(req: DetectTypeRequest) -> Dict[str, Any]:
    """
    Detect reaction type using rxn-insight when available, with router fallback.
    
    This function:
    1. Normalizes the reaction SMILES
    2. Attempts rxn-insight detection (if available)
    3. Falls back to router-based family detection
    4. Returns combined results with selected family
    
    Args:
        req: DetectTypeRequest with reaction SMILES
        
    Returns:
        Dictionary containing:
        - input: Normalized reaction SMILES
        - rxn_insight_available: Whether rxn-insight is installed
        - rxn_insight: rxn-insight detection results (or None)
        - router_fallback: Router-based family detection
        - selected_family: Best family determination
        
    Raises:
        ValidationError: If reaction SMILES is invalid
    """
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")
    
    rxn = req.reaction
    
    # Normalize reaction SMILES
    norm = chem.smiles.normalize_reaction(rxn)
    
    # Extract reactants for router fallback
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]
    
    # Router-based family detection (always available)
    fallback = chem.router.detect_family(reactants)
    
    # Try rxn-insight if available
    auto = None
    if _HAS_RXN_INSIGHT:
        try:
            auto = rxn_detect_type(norm.get("normalized") or rxn)
        except Exception:
            # Silently fall back to router if rxn-insight fails
            auto = None
    
    # Select best family
    selected = None
    if isinstance(auto, dict) and (auto.get("mapped_family") or auto.get("success")):
        selected = auto.get("mapped_family") or fallback.get("family")
    else:
        selected = fallback.get("family")
    
    return {
        "input": {"reaction_smiles": norm.get("normalized") or rxn},
        "rxn_insight_available": bool(_HAS_RXN_INSIGHT),
        "rxn_insight": auto,
        "router_fallback": fallback,
        "selected_family": selected,
    }


def is_rxn_insight_available() -> bool:
    """Check if rxn-insight is available for advanced reaction type detection."""
    return _HAS_RXN_INSIGHT
