"""
Normalization helpers for chemical entries, conditions, and recommendations.

This module provides internal helper functions to normalize various data structures
into the standard output schema format.
"""

from __future__ import annotations

import copy
import re
from typing import Any, Dict, List, Optional, Tuple
from chemtools import reagent


_REQUIRED_CHEM_FIELDS: Tuple[str, ...] = (
    "name",
    "abbreviation",
    "cas",
    "smiles",
    "equivalents",
    "role",
)


def normalize_chemical_entry(entry: Dict[str, Any], role_hint: Optional[str] = None) -> Dict[str, Any]:
    """
    Ensure chemical entries provide the required schema keys.
    
    Args:
        entry: Chemical entry dictionary
        role_hint: Optional role to use if not present in entry
        
    Returns:
        Normalized chemical dictionary with all required fields
    """
    data: Dict[str, Any] = copy.deepcopy(entry or {})
    
    if role_hint and not data.get("role"):
        data["role"] = role_hint
    
    for key in _REQUIRED_CHEM_FIELDS:
        data.setdefault(key, None)
    
    return data


def normalize_condition_value(value: Any, default_unit: Optional[str] = None) -> Any:
    """
    Normalize temperature/time values to a consistent structure.
    
    Handles various input formats:
    - Numbers: converted to {value: X, unit: default_unit}
    - Dicts: normalized with 'value' and 'unit' keys
    - Strings: parsed for numeric value and unit
    
    Args:
        value: Condition value (number, string, or dict)
        default_unit: Default unit to use if not specified
        
    Returns:
        Normalized condition value or None
    """
    if value is None:
        return None
    
    if isinstance(value, dict):
        normalized = copy.deepcopy(value)
        if "value" not in normalized:
            if "amount" in normalized:
                normalized["value"] = normalized.pop("amount")
            elif "quantity" in normalized:
                normalized["value"] = normalized.pop("quantity")
        if default_unit and "unit" not in normalized and "units" not in normalized:
            normalized.setdefault("unit", default_unit)
        return normalized
    
    if isinstance(value, (int, float)):
        result: Dict[str, Any] = {"value": value}
        if default_unit:
            result["unit"] = default_unit
        return result
    
    if isinstance(value, str):
        match = re.match(r"\s*([+-]?\d+(?:\.\d+)?)", value)
        if match:
            numeric = float(match.group(1))
            remainder = value[match.end():].strip()
            result = {"value": numeric}
            unit = remainder or default_unit
            if unit:
                result["unit"] = unit
            result["text"] = value
            return result
        return value
    
    return value


def normalize_conditions_block(conditions: Any) -> Dict[str, Any]:
    """
    Normalize conditions dictionaries into the standard schema.
    
    Handles temperature, time, and atmosphere, plus any extra metadata.
    
    Args:
        conditions: Conditions dictionary or value
        
    Returns:
        Normalized conditions dictionary with standard keys
    """
    if not isinstance(conditions, dict):
        return {
            "temperature": None,
            "time": None,
            "atmosphere": None,
            "note": conditions,
        }
    
    normalized: Dict[str, Any] = {}
    
    # Temperature
    temp_value = (
        conditions.get("temperature")
        or conditions.get("temperature_C")
        or conditions.get("temp_C")
        or conditions.get("temp")
    )
    normalized["temperature"] = normalize_condition_value(temp_value, "C")
    
    # Time
    time_value = (
        conditions.get("time")
        or conditions.get("time_h")
        or conditions.get("time_hr")
        or conditions.get("time_hours")
    )
    normalized["time"] = normalize_condition_value(time_value, "h")
    
    # Atmosphere
    normalized["atmosphere"] = conditions.get("atmosphere") or conditions.get("environment")
    
    # Preserve any additional condition metadata
    extras = {
        key: copy.deepcopy(value)
        for key, value in conditions.items()
        if key not in {"temperature", "temperature_C", "temp_C", "temp", "time", "time_h", "time_hr", "time_hours", "atmosphere", "environment"}
    }
    if extras:
        normalized["extras"] = extras
    
    for key in ("temperature", "time", "atmosphere"):
        normalized.setdefault(key, None)
    
    return normalized


def normalize_recommendation_entry(rec: Dict[str, Any], position: int) -> Dict[str, Any]:
    """
    Normalize a single recommendation entry to the standard schema.
    
    Args:
        rec: Recommendation entry dictionary
        position: Position/rank in the recommendations list
        
    Returns:
        Normalized recommendation dictionary
    """
    source = rec or {}
    normalized: Dict[str, Any] = {}
    
    normalized["rank"] = source.get("rank") or position
    
    chemicals_source = source.get("chemicals")
    if not isinstance(chemicals_source, list) or not chemicals_source:
        chemicals_source = source.get("reagents") or []
    
    normalized_chemicals: List[Dict[str, Any]] = []
    for chem in chemicals_source:
        if isinstance(chem, dict):
            normalized_chemicals.append(normalize_chemical_entry(chem))
        else:
            normalized_chemicals.append(normalize_chemical_entry({"name": str(chem)}, None))
    normalized["chemicals"] = normalized_chemicals
    
    normalized["conditions"] = normalize_conditions_block(source.get("conditions", {}))
    
    for key in ("confidence", "precedent_count", "summary", "component_scores", "reasons", "source", "reaction"):
        if key in source:
            normalized[key] = copy.deepcopy(source[key])
    
    extras_keys = set(source.keys()) - {
        "rank",
        "chemicals",
        "reagents",
        "conditions",
        "confidence",
        "precedent_count",
        "summary",
        "component_scores",
        "reasons",
        "source",
        "reaction",
    }
    if extras_keys:
        normalized["extras"] = {key: copy.deepcopy(source[key]) for key in extras_keys}
    
    return normalized


def normalize_recommendations(recommendations_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Normalize a list of recommendation entries.
    
    Args:
        recommendations_data: List of recommendation dictionaries
        
    Returns:
        List of normalized recommendation dictionaries
    """
    if not recommendations_data:
        return []
    return [
        normalize_recommendation_entry(rec, position=index)
        for index, rec in enumerate(recommendations_data, start=1)
    ]


def parse_amount_to_equivalents(amount: str, role: str) -> Optional[float]:
    """
    Parse amount string to equivalents.
    
    Handles formats like:
    - "200%" or "200–250%" -> 2.0 or 2.25 (midpoint)
    - "0.5–1.5%" -> 0.01 (midpoint)
    - "2.0 eq" -> 2.0
    - "10–20%" for ligands -> 0.15 (midpoint)
    
    Args:
        amount: Amount string from rule file
        role: Reagent role (metal_source, ligand, base, solvent, etc.)
    
    Returns:
        Equivalents as float, or None if cannot parse
    """
    if not amount or not isinstance(amount, str):
        return None
    
    amount = amount.strip()
    
    # Check for explicit "eq" notation
    eq_match = re.search(r'([\d.]+)\s*eq', amount, re.IGNORECASE)
    if eq_match:
        try:
            return float(eq_match.group(1))
        except ValueError:
            return None
    
    # Check for percentage notation
    # Match patterns like "200%", "0.5–1.5%", "10-20%"
    pct_match = re.search(r'([\d.]+)(?:\s*[–\-]\s*([\d.]+))?\s*%', amount)
    if pct_match:
        try:
            low = float(pct_match.group(1))
            high = float(pct_match.group(2)) if pct_match.group(2) else low
            
            # For catalysts/ligands, % usually means mol%
            # For bases, % usually means % equivalents (e.g., 200% = 2.0 eq)
            if role in ['metal_catalyst', 'metal_source', 'catalyst', 'ligand']:
                # Convert mol% to equivalents (e.g., 1.5% = 0.015 eq)
                midpoint = (low + high) / 2.0
                return midpoint / 100.0
            else:
                # For bases and other reagents, % is equivalents (e.g., 200% = 2.0 eq)
                midpoint = (low + high) / 2.0
                return midpoint / 100.0
        except (ValueError, TypeError):
            return None
    
    # Check for "M" notation (molarity - not equivalents)
    if 'M' in amount or 'm' in amount.lower():
        return None
    
    return None


def normalize_rule_string_value(value: str, role: str, amount: Optional[str] = None) -> Dict[str, Any]:
    """
    Convert string-based rule entries into structured chemical dictionaries with database enrichment.
    
    Args:
        value: Reagent name/description from rule file
        role: Chemical role (metal_catalyst, ligand, base, etc.)
        amount: Optional amount string to parse for equivalents
        
    Returns:
        Normalized chemical dictionary with enriched data
    """
    # Remove parenthetical notes to get clean name for lookup
    cleaned = re.sub(r"\s*\(.*?\)\s*$", "", value or "").strip()
    
    # Map output roles to database reagent types
    role_to_type_map = {
        "metal_catalyst": "metal_catalyst",
        "catalyst": "metal_catalyst",
        "ligand": "ligand",
        "base": "base",
        "solvent": "solvent",
        "additive": "additive",
        "oxidant": "oxidant",
        "reductant": "reductant",
        "acid": "acid",
    }
    
    reagent_type = role_to_type_map.get(role, role)
    
    # Try to enrich with database information
    enriched_info = reagent.enrich_reagent_info(cleaned, reagent_type)
    
    # Build chemical entry
    chemical: Dict[str, Any] = {
        "name": enriched_info.get("name") if enriched_info.get("found") else cleaned or value or None,
        "role": role,
        "notes": value,
        "abbreviation": enriched_info.get("abbreviation") if enriched_info.get("found") else None,
        "cas": enriched_info.get("cas") if enriched_info.get("found") else None,
        "smiles": enriched_info.get("smiles") if enriched_info.get("found") else None,
    }
    
    # Parse equivalents from amount parameter (new standardized format)
    equivalents = None
    if amount:
        equivalents = parse_amount_to_equivalents(amount, role)
    
    # Fallback: extract equivalents from the value string itself (legacy format)
    if equivalents is None:
        eq_match = re.search(r"([+-]?\d+(?:\.\d+)?)\s*eq", value or "", re.IGNORECASE)
        if eq_match:
            try:
                equivalents = float(eq_match.group(1))
            except ValueError:
                pass
    
    # Set equivalents if found
    if equivalents is not None:
        chemical["equivalents"] = equivalents
    else:
        chemical["equivalents"] = None
    
    # Extract mol% from the original value string for additional info
    mol_match = re.search(r"([+-]?\d+(?:\.\d+)?)\s*mol%?", value or "", re.IGNORECASE)
    if mol_match:
        try:
            chemical.setdefault("extras", {})["mol_percent"] = float(mol_match.group(1))
        except ValueError:
            pass
    
    return normalize_chemical_entry(chemical, role)
