"""
Utility functions for output formatting.

This module contains miscellaneous helper functions for enriching reagent
information, parsing condition strings, and other formatting utilities.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Optional, Tuple
from chemtools import reagent


def enrich_reagent(
    name: str,
    reagent_type: str,
    role: str,
    equivalents: Optional[float] = None,
    equiv_range: Optional[Tuple[float, float]] = None,
    smiles: Optional[str] = None,
    reagent_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Enrich a single reagent with database information.
    
    Args:
        name: Reagent name or abbreviation
        reagent_type: Type for database lookup (ligand, base, solvent, etc.)
        role: Reagent role (catalyst, ligand, base, solvent, etc.)
        equivalents: Stoichiometric equivalents
        equiv_range: Range of equivalents [min, max]
        smiles: SMILES string (for starting materials)
        reagent_id: Optional unique identifier
    
    Returns:
        Enriched reagent dictionary
    """
    # Try to enrich from database
    enriched = {
        "id": reagent_id or f"{role.upper()[:3]}_1",
        "role": role,
    }
    
    # For starting materials, use SMILES
    if role in ["electrophile", "nucleophile", "starting_material"] and smiles:
        enriched.update({
            "name": None,
            "abbreviation": None,
            "cas": None,
            "smiles": smiles,
            "inchi_key": None,
        })
    else:
        # Look up in database
        db_info = reagent.enrich_reagent_info(name, reagent_type)
        
        if db_info.get("found"):
            enriched.update({
                "name": db_info.get("name"),
                "abbreviation": db_info.get("abbreviation"),
                "cas": db_info.get("cas"),
                "smiles": db_info.get("smiles"),
                "inchi_key": db_info.get("inchi_key"),
            })
            
            # Add category for catalysts
            if role == "catalyst":
                enriched["category"] = reagent_type
            
            # Add properties from database
            if db_info.get("roles") and reagent_type in db_info["roles"]:
                role_props = db_info["roles"][reagent_type]
                enriched["properties"] = role_props
        else:
            # Not found in database - use provided info
            enriched.update({
                "name": name,
                "abbreviation": None,
                "cas": None,
                "smiles": smiles,
                "inchi_key": None,
            })
    
    # Add equivalents
    if equivalents is not None:
        if equiv_range:
            enriched["equivalents"] = {
                "value": equivalents,
                "range": list(equiv_range),
                "unit": "eq"
            }
        else:
            enriched["equivalents"] = {
                "value": equivalents,
                "range": [equivalents, equivalents],
                "unit": "eq"
            }
    
    # Add loading for catalysts and ligands
    if role in ["catalyst", "ligand"] and equivalents is not None:
        loading_pct = equivalents * 100
        enriched["loading"] = {
            "value": round(loading_pct, 1),
            "range": [round(equiv_range[0] * 100, 1), round(equiv_range[1] * 100, 1)] if equiv_range else [loading_pct, loading_pct],
            "unit": "mol%"
        }
    
    return enriched


def format_conditions(
    temperature: Optional[float] = None,
    temp_range: Optional[Tuple[float, float]] = None,
    time_hours: Optional[float] = None,
    time_range: Optional[Tuple[float, float]] = None,
    atmosphere: Optional[str] = None,
    pressure: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Format reaction conditions.
    
    Args:
        temperature: Temperature value (°C)
        temp_range: Temperature range [min, max]
        time_hours: Reaction time (hours)
        time_range: Time range [min, max]
        atmosphere: Atmosphere description
        pressure: Pressure (atm)
        
    Returns:
        Formatted conditions dictionary
    """
    conditions = {}
    
    if temperature is not None:
        conditions["temperature"] = {
            "value": temperature,
            "range": list(temp_range) if temp_range else [temperature, temperature],
            "unit": "°C"
        }
    
    if time_hours is not None:
        conditions["time"] = {
            "value": time_hours,
            "range": list(time_range) if time_range else [time_hours, time_hours],
            "unit": "hours"
        }
    
    if atmosphere:
        conditions["atmosphere"] = {
            "type": "inert",
            "gas": atmosphere,
            "requirement": "anhydrous" if "N" in atmosphere or "Ar" in atmosphere else "ambient"
        }
    
    if pressure:
        conditions["pressure"] = {
            "value": pressure,
            "unit": "atm"
        }
    
    return conditions


def format_recommendation(
    rank: int,
    confidence: float,
    support: int,
    reaction_smiles: str,
    reagents: list,
    conditions: dict,
    similarity_score: Optional[float] = None,
    expected_yield: Optional[float] = None,
    yield_range: Optional[Tuple[float, float]] = None,
    precedents: Optional[list] = None,
    reasoning: Optional[dict] = None,
    notes: Optional[str] = None,
    tags: Optional[list] = None,
) -> Dict[str, Any]:
    """Format a single recommendation with enhanced metadata."""
    recommendation = {
        "rank": rank,
        "confidence": round(confidence, 4),
        "support": support,
        "reaction": {"smiles": reaction_smiles},
        "reagents": reagents,
        "conditions": conditions,
    }
    
    if similarity_score is not None:
        recommendation["similarity_score"] = round(similarity_score, 4)
    
    if expected_yield is not None:
        recommendation["expected_outcome"] = {
            "yield": {
                "typical": round(expected_yield, 1),
                "range": [round(yield_range[0], 1), round(yield_range[1], 1)] if yield_range else [expected_yield, expected_yield],
                "unit": "%"
            }
        }
    
    if precedents:
        recommendation["precedents"] = {
            "count": support,
            "top_similar": precedents[:5]
        }
    
    if reasoning:
        recommendation["reasoning"] = reasoning
    elif support > 0:
        recommendation["reasoning"] = {
            "method": "precedent-based",
            "basis": f"Based on {support} similar precedent reactions",
            "confidence_factors": [
                f"High structural similarity (score: {similarity_score:.3f})" if similarity_score and similarity_score > 0.8 else None,
                f"{support} successful precedent reactions" if support > 10 else None,
                "Standard conditions for this reaction type",
            ],
        }
        recommendation["reasoning"]["confidence_factors"] = [
            f for f in recommendation["reasoning"]["confidence_factors"] if f is not None
        ]
    
    if notes:
        recommendation["notes"] = notes
    
    if tags:
        recommendation["tags"] = tags
    
    return recommendation


def parse_condition_string(condition_str: str) -> Tuple[str, Optional[float]]:
    """
    Parse condition strings like 'K2CO3 2.0 eq' into name and equivalents.
    
    Args:
        condition_str: Condition string
    
    Returns:
        Tuple of (name, equivalents)
    """
    # Try to extract equivalents
    match = re.search(r'(\d+\.?\d*)\s*eq', condition_str, re.IGNORECASE)
    
    if match:
        equiv = float(match.group(1))
        # Remove equiv part to get name
        name = re.sub(r'\s*\d+\.?\d*\s*eq.*$', '', condition_str, flags=re.IGNORECASE).strip()
        return name, equiv
    
    return condition_str.strip(), None
