"""
Simplified core output formatters for molecule and reaction featurization.

This module provides streamlined outputs with 6-8 essential fields for common use cases.
For detailed analysis, use the extended format via detailed=True option.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional


def build_core_molecule(full_molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract core molecule features from full featurization output.
    
    Core output contains:
    - smiles: Molecule SMILES string
    - motifs: Primary functional groups (simplified)
    - properties: Key steric/electronic summary
    - rdkit: Standard molecular descriptors
    
    Args:
        full_molecule: Complete molecule bundle from molecule.py
        
    Returns:
        Simplified core molecule with 4-6 fields
    """
    # Extract simplified motifs (just id and rank)
    motifs = []
    for m in full_molecule.get("motifs", []):
        motif = {
            "id": m.get("compound_id", "unknown"),
            "rank": m.get("rank_score", 0),
        }
        # Include atom indices if available
        if m.get("a_atom_idx") is not None:
            motif["a_atom_idx"] = m.get("a_atom_idx")
        if m.get("b_atom_idx") is not None:
            motif["b_atom_idx"] = m.get("b_atom_idx")
        motifs.append(motif)
    
    # Extract key properties from analyses
    max_steric = 0.0
    electronic_scores = []
    for analysis in full_molecule.get("analyses", []):
        if isinstance(analysis, dict):
            steric = analysis.get("steric", {})
            if isinstance(steric, dict):
                score = steric.get("score_0_10", 0)
                if isinstance(score, (int, float)):
                    max_steric = max(max_steric, float(score))
            
            electronic = analysis.get("electronic")
            if isinstance(electronic, dict):
                score = electronic.get("score_0_10")
                if isinstance(score, (int, float)):
                    electronic_scores.append(float(score))
            elif isinstance(electronic, list):
                for e in electronic:
                    if isinstance(e, dict):
                        score = e.get("score_0_10")
                        if isinstance(score, (int, float)):
                            electronic_scores.append(float(score))
    
    avg_electronic = round(sum(electronic_scores) / len(electronic_scores), 1) if electronic_scores else 5.0
    
    properties = {
        "max_steric": round(max_steric, 1),
        "avg_electronic": avg_electronic,
    }
    
    # Extract RDKit properties if available
    rdkit = full_molecule.get("rdkit_properties", {})
    if rdkit:
        rdkit_core = {
            "mw": round(rdkit.get("molecular_weight", 0), 1),
            "logP": round(rdkit.get("logP", 0), 1),
            "tpsa": round(rdkit.get("TPSA", 0), 1),
        }
    else:
        rdkit_core = {}
    
    core = {
        "smiles": full_molecule.get("smiles", ""),
        "motifs": motifs,
        "properties": properties,
    }
    
    if rdkit_core:
        core["rdkit"] = rdkit_core
    
    return core


def build_core_reaction(full_reaction: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract core reaction features from full featurization output.
    
    Core output contains:
    - reaction_smiles: Reaction SMILES string
    - reaction_type: Best match only (name + confidence)
    - reaction_key: Formatted motif summary
    - reactants: Simplified molecule bundles
    - products: Simplified molecule bundles
    - feasibility: Simple status (high/medium/low/unknown)
    
    Args:
        full_reaction: Complete reaction bundle from reaction.py
        
    Returns:
        Simplified core reaction with 6-7 fields
    """
    # Extract reaction type (best match only)
    reaction_type_full = full_reaction.get("reaction_type", {})
    reaction_type = reaction_type_full.get("name") or reaction_type_full.get("reaction_type") or "Unknown"
    confidence = reaction_type_full.get("confidence", 0.0)
    
    # Simplify reactants and products
    reactants = []
    for r in full_reaction.get("reactants", []):
        reactants.append(build_core_molecule(r))
    
    products = []
    for p in full_reaction.get("products", []):
        products.append(build_core_molecule(p))
    
    # Simplify feasibility to enum
    feasibility_full = full_reaction.get("feasibility", {})
    if feasibility_full:
        is_feasible = feasibility_full.get("feasible", False)
        feas_confidence = feasibility_full.get("confidence", 0.0)
        if is_feasible:
            if feas_confidence >= 0.7:
                feasibility = "high"
            else:
                feasibility = "medium"
        else:
            feasibility = "low"
    else:
        feasibility = "unknown"
    
    core = {
        "reaction_smiles": full_reaction.get("reaction_smiles", ""),
        "reaction_type": reaction_type,
        "confidence": round(confidence, 3),
        "reaction_key": full_reaction.get("reaction_key", ""),
        "reactants": reactants,
        "products": products,
        "feasibility": feasibility,
    }
    
    return core


def build_extended_molecule(full_molecule: Dict[str, Any]) -> Dict[str, Any]:
    """
    Build extended molecule output with detailed analysis.
    
    Extended output includes core fields plus:
    - extended.per_motif_analysis: Detailed steric/electronic/nearby groups
    - extended.snar_feasibility: SNAr reaction feasibility check
    - extended.context_motifs: Background/scaffold motifs
    - extended.ranked_motifs: Sorted by importance
    
    Args:
        full_molecule: Complete molecule bundle from molecule.py
        
    Returns:
        Core fields + extended analysis section
    """
    core = build_core_molecule(full_molecule)
    
    # Build per-motif analysis from full analyses
    per_motif_analysis = []
    for analysis in full_molecule.get("analyses", []):
        if not isinstance(analysis, dict):
            continue
        
        motif_analysis: Dict[str, Any] = {
            "motif_id": analysis.get("compound_id", "unknown"),
        }
        
        # Add steric if available
        steric = analysis.get("steric")
        if isinstance(steric, dict):
            motif_analysis["steric"] = {
                "score": steric.get("score_0_10", 0),
                "classification": steric.get("classification", ""),
                "description": steric.get("description", ""),
            }
        
        # Add electronic if available
        electronic = analysis.get("electronic")
        if isinstance(electronic, dict):
            motif_analysis["electronic"] = {
                "score": electronic.get("score_0_10", 5),
                "description": electronic.get("description", ""),
            }
        elif isinstance(electronic, list) and electronic:
            # Take first electronic analysis
            e = electronic[0]
            if isinstance(e, dict):
                motif_analysis["electronic"] = {
                    "score": e.get("score_0_10", 5),
                    "description": e.get("description", ""),
                }
        
        # Add nearby groups if available
        nearby = analysis.get("nearby_groups", [])
        if nearby:
            motif_analysis["nearby_groups"] = [
                g.get("name") or g.get("group_id", "unknown")
                for g in nearby
                if isinstance(g, dict)
            ]
        
        per_motif_analysis.append(motif_analysis)
    
    extended: Dict[str, Any] = {
        "per_motif_analysis": per_motif_analysis,
    }
    
    # Add SNAr feasibility if available
    snar = full_molecule.get("snar_feasibility", [])
    if snar:
        extended["snar_feasibility"] = snar
    
    # Add context and ranked motifs
    context = full_molecule.get("context_motifs", [])
    if context:
        extended["context_motifs"] = [
            m.get("compound_id", "unknown") for m in context if isinstance(m, dict)
        ]
    
    ranked = full_molecule.get("ranked_motifs", [])
    if ranked:
        extended["ranked_motifs"] = ranked
    
    core["extended"] = extended
    return core


def build_extended_reaction(full_reaction: Dict[str, Any]) -> Dict[str, Any]:
    """
    Build extended reaction output with detailed analysis.
    
    Extended output includes core fields plus:
    - extended.detection: All detection matches (not just best)
    - extended.aggregates: Reaction-wide statistics
    - extended.role_classification: Reactant/agent role analysis
    - extended.normalization_log: Input normalization details
    
    Args:
        full_reaction: Complete reaction bundle from reaction.py
        
    Returns:
        Core fields + extended analysis section
    """
    core = build_core_reaction(full_reaction)
    
    # Replace simplified reactants/products with extended versions
    reactants_extended = []
    for r in full_reaction.get("reactants", []):
        reactants_extended.append(build_extended_molecule(r))
    
    products_extended = []
    for p in full_reaction.get("products", []):
        products_extended.append(build_extended_molecule(p))
    
    core["reactants"] = reactants_extended
    core["products"] = products_extended
    
    # Build extended section
    extended: Dict[str, Any] = {}
    
    # Add detection matches
    detection = full_reaction.get("detection", {})
    if detection:
        matches = detection.get("matches", [])
        if matches:
            # Format matches with top 5
            extended["detection"] = {
                "matches": matches[:5],
                "total_matches": len(matches),
            }
        # Include validation metadata if present
        validation = detection.get("validation")
        if validation:
            if "detection" not in extended:
                extended["detection"] = {}
            extended["detection"]["validation"] = validation
    
    # Add aggregates
    aggregates = full_reaction.get("aggregates", {})
    if aggregates:
        extended["aggregates"] = aggregates
    
    # Add role classification
    roles = full_reaction.get("roles")
    if roles:
        extended["role_classification"] = {
            "reactants": roles,
        }
    
    agent_roles = full_reaction.get("agent_roles")
    if agent_roles:
        if "role_classification" not in extended:
            extended["role_classification"] = {}
        extended["role_classification"]["agents"] = agent_roles
    
    # Add intramolecular info
    intramolecular = full_reaction.get("intramolecular")
    if intramolecular:
        extended["intramolecular"] = intramolecular
    
    # Add normalization log
    normalized = full_reaction.get("normalized", {})
    if normalized:
        extended["normalization_log"] = {
            "input": normalized.get("input"),
            "normalized": normalized.get("normalized"),
            "errors": normalized.get("errors", []),
        }
    
    if extended:
        core["extended"] = extended
    
    return core


__all__ = [
    "build_core_molecule",
    "build_core_reaction",
    "build_extended_molecule",
    "build_extended_reaction",
]
