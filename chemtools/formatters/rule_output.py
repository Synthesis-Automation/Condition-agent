"""
Rule-based output formatting for SchemeConditionDB and pattern-matching systems.

This module handles formatting of rule-based recommendation outputs, including
conversion of SCDB match results into standardized recommendations.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, List, Optional
from .normalization import (
    normalize_chemical_entry,
    normalize_conditions_block,
    normalize_rule_string_value,
)


def starting_material_entries(reaction_smiles: str) -> List[Dict[str, Any]]:
    """
    Build chemical entries for starting materials based on reaction SMILES.
    
    Args:
        reaction_smiles: Reaction SMILES string
        
    Returns:
        List of starting material chemical entries
    """
    if ">>" not in reaction_smiles:
        return []
    
    reactants = [frag for frag in reaction_smiles.split(">>", 1)[0].split(".") if frag]
    roles = ["electrophile", "nucleophile"]
    
    entries = []
    for idx, smiles in enumerate(reactants):
        role = roles[idx] if idx < len(roles) else "starting_material"
        entries.append(normalize_chemical_entry({"smiles": smiles, "role": role}, role))
    
    return entries


def convert_rule_match_to_recommendations(
    reaction_smiles: str,
    match_result: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Convert SCDB match results into standardized recommendation entries.
    
    Supports both new standardized format (with 'reagents' array) and legacy format
    (with direct keys like 'pd_source', 'ligand', 'base', 'solvent').
    
    Args:
        reaction_smiles: Input reaction SMILES
        match_result: Raw SCDB match result dictionary
        
    Returns:
        List of standardized recommendation entries
    """
    conditions = match_result.get("conditions")
    if isinstance(conditions, dict):
        condition_entries = [conditions]
    elif isinstance(conditions, list):
        condition_entries = [cond for cond in conditions if isinstance(cond, dict)]
    else:
        condition_entries = []
    
    if not condition_entries:
        return []
    
    # Check if entry contains reagents data (new format is at match_result level via raw entry)
    entry_raw = match_result.get("trace", {}).get("selected", {})
    
    starting_materials = starting_material_entries(reaction_smiles)
    recommendations: List[Dict[str, Any]] = []
    
    for idx, cond in enumerate(condition_entries, start=1):
        chemicals: List[Dict[str, Any]] = copy.deepcopy(starting_materials)
        
        # Role mapping for new standardized format
        role_mapping = {
            "metal_source": "metal_catalyst",
            "metal_catalyst": "metal_catalyst",
            "catalyst": "metal_catalyst",
            "ligand": "ligand",
            "base": "base",
            "solvent": "solvent",
            "additive": "additive",
            "partner": "partner",
            "boron_partner": "partner",
        }
        
        # Try new format first: check for reagents in the entry's raw data
        reagents_added = False
        if "reagents" in match_result:
            # Reagents at match result level (passed through from entry)
            reagents = match_result.get("reagents", [])
            if isinstance(reagents, list):
                for reagent in reagents:
                    if isinstance(reagent, dict):
                        name = reagent.get("name")
                        role = reagent.get("role", "reagent")
                        amount = reagent.get("amount")
                        
                        # Map role to output role
                        output_role = role_mapping.get(role, role)
                        
                        if name:
                            chemicals.append(normalize_rule_string_value(str(name), output_role, amount=amount))
                            reagents_added = True
        
        # Fallback to legacy format if reagents not found in new format
        if not reagents_added:
            def add_entry(value: Any, role: str) -> None:
                if not value:
                    return
                if isinstance(value, dict):
                    chemicals.append(normalize_chemical_entry(value, role))
                else:
                    chemicals.append(normalize_rule_string_value(str(value), role))
            
            # Legacy format extraction
            add_entry(cond.get("pd_source") or cond.get("catalyst") or cond.get("metal_source"), "metal_catalyst")
            add_entry(cond.get("ligand"), "ligand")
            add_entry(cond.get("base"), "base")
            add_entry(cond.get("solvent"), "solvent")
            add_entry(cond.get("additive"), "additive")
            add_entry(cond.get("boron_partner") or cond.get("partner"), "partner")
            add_entry(cond.get("condition_core"), "condition_core")
        
        recommendation = {
            "rank": idx,
            "chemicals": chemicals,
            "conditions": normalize_conditions_block(cond),
        }
        
        loadings = cond.get("loadings")
        if isinstance(loadings, dict) and loadings:
            recommendation.setdefault("conditions", {}).setdefault("extras", {})["loadings"] = loadings
        
        source_meta = {
            "match_type": match_result.get("match_type"),
            "entry_id": match_result.get("entry_id"),
            "entry_name": match_result.get("entry_name"),
            "priority": match_result.get("priority"),
        }
        if any(value is not None for value in source_meta.values()):
            recommendation["source"] = source_meta
        
        recommendations.append(recommendation)
    
    return recommendations
