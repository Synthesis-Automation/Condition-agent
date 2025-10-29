"""
Matching Service - SMILES normalization and reaction family detection.

This service handles:
- SMILES normalization (canonicalization)
- Reaction family detection using unified detection API
- Reaction type detection using ML-enhanced detection
"""

from typing import Dict, Any, List
from chemtools import chem
from chemtools import detect_reaction  # New unified API
from chemtools.recommend.utils import friendly_family_label
from chemtools.contracts import NormalizeRequest, DetectFamilyRequest, DetectTypeRequest
from chemtools.exceptions import ValidationError


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
    
    Uses the new unified detection API (chemtools.detect_reaction).
    
    Args:
        req: DetectFamilyRequest with reactants list
        
    Returns:
        Dictionary with family detection results (old schema for compatibility)
        
    Raises:
        ValidationError: If reactants are invalid
    """
    if not req.reactants:
        raise ValidationError("Reactants list cannot be empty")
    
    # Convert reactants to pseudo-reaction for unified API
    reaction = ".".join(req.reactants) + ">>"
    
    # Use new unified API with rule-based detection only
    result = detect_reaction(reaction, use_ml=False)
    
    # Convert to old schema for backwards compatibility
    family = result["family"]
    confidence = result["confidence"]
    hits = result["details"].get("functional_groups", {})
    
    response = {
        "family": family,
        "confidence": confidence,
        "hits": hits
    }
    
    # Add friendly label
    label = friendly_family_label(family)
    if label:
        response["family_label"] = label
    
    return response


def detect_reaction_type(req: DetectTypeRequest) -> Dict[str, Any]:
    """
    Detect reaction type using unified detection API with ML enhancement.
    
    This function uses the new detect_reaction() API which combines:
    1. Rule-based detection (SMARTS patterns)
    2. ML detection (rxn-insight if available)
    3. Catalyst-aware refinements
    
    Args:
        req: DetectTypeRequest with reaction SMILES
        
    Returns:
        Dictionary containing:
        - input: Normalized reaction SMILES
        - detection: Full detection results from unified API
        - selected_family: Best family determination
        - ml_available: Whether ML detection was used
        
    Raises:
        ValidationError: If reaction SMILES is invalid
    """
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")
    
    rxn = req.reaction
    
    # Normalize reaction SMILES
    norm = chem.smiles.normalize_reaction(rxn)
    normalized_rxn = norm.get("normalized") or rxn
    
    # Use new unified detection API with ML enabled
    detection_result = detect_reaction(normalized_rxn, use_ml=True)
    
    # Extract key fields
    family = detection_result["family"]
    confidence = detection_result["confidence"]
    method = detection_result["method"]
    details = detection_result["details"]
    
    # Check if ML was actually used
    ml_used = "ml" in method.lower() or details.get("ml_prediction") is not None
    
    # Add friendly labels
    family_label = friendly_family_label(family)
    
    # Build response with new and legacy fields
    response = {
        "input": {"reaction_smiles": normalized_rxn},
        "detection": {
            "family": family,
            "confidence": confidence,
            "method": method,
            "family_label": family_label,
        },
        "selected_family": family,
        "selected_family_label": family_label,
        "ml_available": ml_used,
        "details": details,  # Full details for debugging
    }
    
    # Add ML prediction details if available
    if details.get("ml_prediction"):
        ml_pred = details["ml_prediction"]
        response["ml_prediction"] = {
            "rxn_class": ml_pred.get("rxn_class"),
            "rxn_name": ml_pred.get("rxn_name"),
            "confidence": ml_pred.get("confidence"),
            "available": ml_pred.get("available", True),
        }
    
    # Add rule-based prediction for comparison
    if details.get("rule_prediction"):
        rule_pred = details["rule_prediction"]
        response["rule_prediction"] = {
            "family": rule_pred.get("family"),
            "confidence": rule_pred.get("confidence"),
        }
    
    return response


def is_rxn_insight_available() -> bool:
    """
    Check if ML detection is available.
    
    Note: This is a legacy function. The new detect_reaction() API
    automatically handles ML availability internally.
    """
    try:
        from chemtools._ml_helpers import is_available
        return is_available()
    except ImportError:
        return False

