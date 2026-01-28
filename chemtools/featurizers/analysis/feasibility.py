"""
Feasibility analysis for specific reaction types.
All feasibility criteria (compound types, thresholds, etc.) are
defined in taxonomy (chemtools/taxonomy/data/reaction_types.v*.json),
not hardcoded."""

from __future__ import annotations
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple


@lru_cache(maxsize=1)
def _load_snar_config() -> Tuple[List[str], Dict[str, float], Dict[str, str]]:
    """Load SNAr feasibility configuration from taxonomy.
    
    Returns:
        Tuple of (electrophile_types, thresholds, interpretations)
    """
    from ...taxonomy import reaction_catalog
    
    definitions, _ = reaction_catalog.load_reaction_catalog()
    snar_def = definitions.get('Snar_cn')
    
    if not snar_def or not snar_def.metadata:
        # Fallback to safe defaults if taxonomy missing
        return (
            ["Ar-Cl", "Ar-Br", "Ar-I", "Ar-F", "AromN-Cl", "AromN-Br", "AromN-I", "AromN-F"],
            {"neutral_score": 5.0, "minimum_activation": 6.0, "high_activation": 7.0},
            {}
        )
    
    feasibility = snar_def.metadata.get('feasibility', {})
    electronic = feasibility.get('electronic_activation', {})
    
    electrophile_types = feasibility.get('electrophile_types', [])
    thresholds = electronic.get('thresholds', {})
    interpretations = electronic.get('interpretation', {})
    
    return electrophile_types, thresholds, interpretations


def analyze_snar_feasibility(reaction_payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze if an SNAr reaction is likely to be feasible based on electronics.
    
    Uses taxonomy-defined compound types and thresholds from reaction_types.v4.0.json.
    """
    # Load SNAr configuration from taxonomy
    electrophile_types, thresholds, interpretations = _load_snar_config()
    
    reactants = reaction_payload.get("reactants", [])
    
    # Find the electrophile (aryl halide) using taxonomy-defined types
    electrophile = None
    target_motif = None
    for r in reactants:
        motifs = r.get("motifs", [])
        for m in motifs:
            cid = m.get("compound_id", "")
            # Check against taxonomy-defined electrophile types
            if any(cid.startswith(etype) for etype in electrophile_types):
                electrophile = r
                target_motif = m
                break
        if electrophile:
            break
            
    if not electrophile:
        return {
            "feasible": False,
            "reason": "No aryl halide electrophile detected",
            "score": 0.0
        }
        
    # Check electronics of the ipso carbon
    electronics = electrophile.get("electronics", {}).get("aryl", [])
    target_cid = target_motif.get("compound_id")
    
    target_elec = None
    for e in electronics:
        if e.get("compound_id") == target_cid:
            target_elec = e
            break
            
    if not target_elec:
        # Fallback: try to find any electronic score for this reactant
        if electronics:
            target_elec = electronics[0]
            
    if not target_elec:
        return {
            "feasible": False,
            "reason": "Could not determine electronics for the reaction center",
            "score": 0.0
        }
        
    # The 'score_0_10' in electronics is the scaffold score (excluding the halide itself)
    scaffold_score = target_elec.get("result", {}).get("score_0_10", thresholds.get("neutral_score", 5.0))
    
    # Use taxonomy-defined thresholds
    min_activation = thresholds.get("minimum_activation", 6.0)
    high_activation = thresholds.get("high_activation", 7.0)
    neutral_score = thresholds.get("neutral_score", 5.0)
    
    is_feasible = scaffold_score >= min_activation
    confidence = "high" if scaffold_score >= high_activation else "medium" if scaffold_score >= min_activation else "low"
    
    # Build reason using taxonomy-defined interpretations (with fallbacks)
    reason = f"Aryl ring electronics score is {scaffold_score} ({neutral_score}=neutral, >{min_activation}=activated)"
    if scaffold_score < min_activation:
        reason += ". " + interpretations.get("below_6.0", "Ring is likely too electron-rich for SNAr")
    elif scaffold_score >= high_activation:
        reason += ". " + interpretations.get("above_7.0", "Ring is highly activated for SNAr")
    else:
        reason += ". " + interpretations.get("6.0_to_7.0", "Ring is moderately activated for SNAr")
    
    return {
        "feasible": is_feasible,
        "confidence": confidence,
        "score": scaffold_score,
        "reason": reason,
        "details": {
            "scaffold_electronic_score": scaffold_score,
            "target_motif": target_cid or "unknown"
        }
    }


def analyze_molecule_snar_feasibility(molecule_payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Analyze SNAr feasibility for all aryl halides in a molecule.
    """
    # Load SNAr configuration from taxonomy
    electrophile_types, thresholds, interpretations = _load_snar_config()
    
    results = []
    motifs = molecule_payload.get("motifs", [])
    electronics = molecule_payload.get("electronics", {}).get("aryl", [])
    
    for m in motifs:
        cid = m.get("compound_id", "")
        if cid.startswith(("Ar-Cl", "Ar-Br", "Ar-I", "Ar-F", "AromN-Cl", "AromN-Br", "AromN-I", "AromN-F")):
            # Find corresponding electronic score
            target_elec = next((e for e in electronics if e.get("compound_id") == cid), None)
            if not target_elec:
                continue
            
            # Use taxonomy-defined thresholds
            min_activation = thresholds.get("minimum_activation", 6.0)
            high_activation = thresholds.get("high_activation", 7.0)
            neutral_score = thresholds.get("neutral_score", 5.0)
            
            scaffold_score = target_elec.get("result", {}).get("score_0_10", neutral_score)
            is_feasible = scaffold_score >= min_activation
            confidence = "high" if scaffold_score >= high_activation else "medium" if scaffold_score >= min_activation else "low"
            
            # Build reason using taxonomy-defined interpretations
            reason = f"Aryl ring electronics score is {scaffold_score} ({neutral_score}=neutral, >{min_activation}=activated)"
            if scaffold_score < min_activation:
                reason += ". " + interpretations.get("below_6.0", "Ring is likely too electron-rich for SNAr")
            elif scaffold_score >= high_activation:
                reason += ". " + interpretations.get("above_7.0", "Ring is highly activated for SNAr")
            else:
                reason += ". " + interpretations.get("6.0_to_7.0", "Ring is moderately activated for SNAr")
                
            results.append({
                "motif": cid,
                "feasible": is_feasible,
                "confidence": confidence,
                "score": scaffold_score,
                "reason": reason
            })
            
    return results
