"""
Precedent detail building and statistics.

This module provides functions for building detailed precedent summaries
and calculating statistics from precedent data.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional


def build_precedent_details(
    precs: List[Dict[str, Any]],
    chosen_core: Optional[str],
    group: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Build detailed precedent information for the recommendation output.
    
    Args:
        precs: All precedents from the search
        chosen_core: The chosen catalyst core
        group: Precedents matching the chosen core
        
    Returns:
        Dictionary with comprehensive precedent details
    """
    # Get top 10 precedents overall
    top_precedents = []
    for i, p in enumerate(precs[:10], 1):
        precedent_info = {
            "rank": i,
            "reaction_id": p.get("reaction_id"),
            "reaction_smiles": p.get("reaction_smiles"),
            "core": p.get("core") or p.get("condition_core"),
            "yield": p.get("yield"),
            "rxn_type": p.get("rxn_type"),  # Add reaction family/type for dataset identification
        }
        
        # Add detailed reagent information
        catalytic_system = p.get("catalytic_system", [])
        if catalytic_system:
            precedent_info["catalysts"] = [
                {
                    "name": cat.get("name"),
                    "cas": cat.get("cas"),
                    "role": "catalyst/ligand"
                }
                for cat in catalytic_system
            ]
        
        # Add reagents (bases, additives, etc.)
        reagents = p.get("reagents", [])
        if reagents:
            precedent_info["reagents"] = [
                {
                    "name": r.get("name"),
                    "cas": r.get("cas"),
                    "role": r.get("role", "reagent").lower()
                }
                for r in reagents
            ]
        
        # Add solvents
        solvents = p.get("solvents", [])
        if solvents:
            precedent_info["solvents"] = [
                {
                    "name": s.get("name"),
                    "cas": s.get("cas"),
                }
                for s in solvents
            ]
        
        # Add conditions
        conditions_data = {}
        if p.get("T_C") is not None:
            conditions_data["temperature_C"] = p.get("T_C")
        if p.get("time_h") is not None:
            conditions_data["time_h"] = p.get("time_h")
        
        # Add from conditions dict if available
        cond_dict = p.get("conditions", {})
        if isinstance(cond_dict, dict):
            if "temperature_c" in cond_dict:
                conditions_data["temperature_C"] = cond_dict["temperature_c"]
            if "time_h" in cond_dict:
                conditions_data["time_h"] = cond_dict["time_h"]
            if "yield_pct" in cond_dict:
                conditions_data["yield_pct"] = cond_dict["yield_pct"]
        
        if conditions_data:
            precedent_info["conditions"] = conditions_data
        
        # Add reference if available
        if p.get("reference"):
            precedent_info["reference"] = p.get("reference")
        
        top_precedents.append(precedent_info)
    
    # Get precedents matching the chosen core
    core_precedents = []
    if chosen_core and group:
        for i, p in enumerate(group[:5], 1):
            core_prec = {
                "rank": i,
                "reaction_id": p.get("reaction_id"),
                "core": p.get("core") or p.get("condition_core"),
                "base": p.get("base_uid"),
                "solvent": p.get("solvent_uid"),
                "yield": p.get("yield"),
            }
            
            # Add temperature and time if available
            if p.get("T_C") is not None:
                core_prec["temperature_C"] = p.get("T_C")
            if p.get("time_h") is not None:
                core_prec["time_h"] = p.get("time_h")
            
            core_precedents.append(core_prec)
    
    return {
        "total_count": len(precs),
        "top_precedents": top_precedents,
        "core_matched_precedents": {
            "core": chosen_core,
            "count": len(group) if group else 0,
            "examples": core_precedents,
        },
        "statistics": {
            "average_yield": calculate_average_yield(precs),
            "yield_range": calculate_yield_range(precs),
            "temperature_range": calculate_temp_range(precs),
            "time_range": calculate_time_range(precs),
        }
    }


def calculate_average_yield(precs: List[Dict[str, Any]]) -> Optional[float]:
    """Calculate average yield from precedents."""
    yields = [p.get("yield") for p in precs if p.get("yield") is not None and isinstance(p.get("yield"), (int, float))]
    if not yields:
        return None
    return round(sum(yields) / len(yields), 1)


def calculate_yield_range(precs: List[Dict[str, Any]]) -> Optional[List[float]]:
    """Calculate yield range from precedents."""
    yields = [p.get("yield") for p in precs if p.get("yield") is not None and isinstance(p.get("yield"), (int, float))]
    if not yields:
        return None
    return [round(min(yields), 1), round(max(yields), 1)]


def calculate_temp_range(precs: List[Dict[str, Any]]) -> Optional[List[float]]:
    """Calculate temperature range from precedents."""
    temps = [p.get("T_C") for p in precs if p.get("T_C") is not None and isinstance(p.get("T_C"), (int, float))]
    if not temps:
        return None
    return [round(min(temps), 1), round(max(temps), 1)]


def calculate_time_range(precs: List[Dict[str, Any]]) -> Optional[List[float]]:
    """Calculate time range from precedents."""
    times = [p.get("time_h") for p in precs if p.get("time_h") is not None and isinstance(p.get("time_h"), (int, float))]
    if not times:
        return None
    return [round(min(times), 1), round(max(times), 1)]
