"""
Feasibility analysis for specific reaction types.
"""

from __future__ import annotations
from typing import Any, Dict, List, Optional

def analyze_snar_feasibility(reaction_payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze if an SNAr reaction is likely to be feasible based on electronics.
    """
    reactants = reaction_payload.get("reactants", [])
    
    # Find the electrophile (aryl halide)
    electrophile = None
    target_motif = None
    for r in reactants:
        motifs = r.get("motifs", [])
        for m in motifs:
            cid = m.get("compound_id", "")
            # Check for aryl halides
            if cid.startswith(("Ar-Cl", "Ar-Br", "Ar-I", "Ar-F", "AromN-Cl", "AromN-Br", "AromN-I", "AromN-F")):
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
    scaffold_score = target_elec.get("result", {}).get("score_0_10", 5.0)
    
    # Thresholds for SNAr
    # 5.0 is neutral (benzene)
    # > 6.0 is starting to be activated
    # > 7.0 is well activated (e.g. nitro group)
    
    is_feasible = scaffold_score >= 6.0
    confidence = "high" if scaffold_score >= 7.0 else "medium" if scaffold_score >= 6.0 else "low"
    
    reason = f"Aryl ring electronics score is {scaffold_score} (5.0=neutral, >6.0=activated)"
    if scaffold_score < 6.0:
        reason += ". Ring is likely too electron-rich for SNAr without strong activation."
    elif scaffold_score >= 7.0:
        reason += ". Ring is highly activated for SNAr."
    else:
        reason += ". Ring is moderately activated for SNAr."
    
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
                
            scaffold_score = target_elec.get("result", {}).get("score_0_10", 5.0)
            is_feasible = scaffold_score >= 6.0
            confidence = "high" if scaffold_score >= 7.0 else "medium" if scaffold_score >= 6.0 else "low"
            
            reason = f"Aryl ring electronics score is {scaffold_score} (5.0=neutral, >6.0=activated)"
            if scaffold_score < 6.0:
                reason += ". Ring is likely too electron-rich for SNAr."
            elif scaffold_score >= 7.0:
                reason += ". Ring is highly activated for SNAr."
            else:
                reason += ". Ring is moderately activated for SNAr."
                
            results.append({
                "motif": cid,
                "feasible": is_feasible,
                "confidence": confidence,
                "score": scaffold_score,
                "reason": reason
            })
            
    return results
