"""
Reagent database lookup utilities.

Provides functions to enrich condition recommendations with detailed reagent information
from the JSON databases in data/reagent_db/.
"""

import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional
from functools import lru_cache


# Cache loaded reagent databases
_REAGENT_CACHE: Dict[str, List[Dict[str, Any]]] = {}


def get_data_dir() -> Path:
    """Get the data directory path."""
    module_dir = Path(__file__).parent.parent.parent
    return module_dir / "data"


@lru_cache(maxsize=32)
def load_reagent_database(reagent_type: str) -> List[Dict[str, Any]]:
    """
    Load a reagent database JSON file.
    
    Args:
        reagent_type: Type of reagent (e.g., 'ligand', 'base', 'solvent', 'metal_catalyst')
        
    Returns:
        List of reagent dictionaries
    """
    if reagent_type in _REAGENT_CACHE:
        return _REAGENT_CACHE[reagent_type]
    
    reagent_file = get_data_dir() / "reagent_db" / f"{reagent_type}.json"
    
    if not reagent_file.exists():
        return []
    
    try:
        with open(reagent_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
            _REAGENT_CACHE[reagent_type] = data
            return data
    except Exception as e:
        print(f"Warning: Failed to load {reagent_file}: {e}")
        return []


def normalize_name(name: str) -> str:
    """Normalize a reagent name for matching."""
    if not name:
        return ""
    # Lowercase, remove extra whitespace, remove common parentheticals
    normalized = name.lower().strip()
    normalized = normalized.replace("\u00b7", ".")
    normalized = re.sub(r"\([^)]*\)", "", normalized)
    normalized = re.sub(r"\[[^]]*\]", "", normalized)
    normalized = re.sub(r"\b(anhydrous|dry|aqueous|aq|solution|soln)\b", "", normalized)
    normalized = re.sub(r"\b\d+(?:\.\d+)?\s*m\b", "", normalized)
    normalized = re.sub(r"(\.|\u00b7)\s*h2o$", "", normalized)
    normalized = ' '.join(normalized.split())  # Collapse whitespace
    # Remove common formatting differences
    normalized = normalized.replace(',', '').replace('-', '').replace('_', '')
    return normalized


def _split_reagent_name(name: str) -> List[str]:
    if not name:
        return []
    parts = re.split(r"[;/]", name)
    return [part.strip() for part in parts if part and part.strip()]


def find_reagent(name: str, reagent_type: str) -> Optional[Dict[str, Any]]:
    """
    Find a reagent by name in the specified database.
    
    Args:
        name: Reagent name to search for
        reagent_type: Type of reagent database to search
        
    Returns:
        Reagent dictionary if found, None otherwise
    """
    if not name:
        return None
    
    db = load_reagent_database(reagent_type)
    if not db:
        return None
    
    # Normalize search name
    search_name = normalize_name(name)
    
    # Try exact match first
    for reagent in db:
        if normalize_name(reagent.get('name', '')) == search_name:
            return reagent
        
        # Check abbreviations
        for abbr in reagent.get('abbreviation', []):
            if normalize_name(abbr) == search_name:
                return reagent
        
        # Check aliases
        for alias in reagent.get('aliases', []):
            if normalize_name(alias) == search_name:
                return reagent
        
        # Check CAS number
        if reagent.get('cas') and normalize_name(reagent.get('cas', '')) == search_name:
            return reagent
    
    # Try partial match
    for reagent in db:
        reagent_name = normalize_name(reagent.get('name', ''))
        if search_name in reagent_name or reagent_name in search_name:
            return reagent

    # Try splitting composite names (e.g., "DMF/H2O")
    if any(sep in name for sep in ["/", ";"]):
        for part in _split_reagent_name(name):
            if part and part != name:
                hit = find_reagent(part, reagent_type)
                if hit:
                    return hit
    
    return None


def enrich_reagent_info(name: str, reagent_type: str) -> Dict[str, Any]:
    """
    Enrich a reagent name with detailed information from database.
    
    Args:
        name: Reagent name
        reagent_type: Type of reagent (ligand, base, solvent, etc.)
        
    Returns:
        Dictionary with reagent details
    """
    # Default response
    result = {
        "name": name,
        "type": reagent_type,
        "cas": None,
        "abbreviation": None,
        "smiles": None,
        "inchi_key": None,
        "found": False,
    }
    
    reagent = find_reagent(name, reagent_type)
    
    if reagent:
        result.update({
            "name": reagent.get('name', name),
            "cas": reagent.get('cas'),
            "abbreviation": reagent.get('abbreviation', [None])[0] if reagent.get('abbreviation') else None,
            "smiles": reagent.get('smiles'),
            "inchi_key": reagent.get('inchi_key'),
            "aliases": reagent.get('aliases', []),
            "roles": reagent.get('roles', {}),
            "found": True,
        })
    
    return result


def enrich_conditions(conditions: Dict[str, Any]) -> Dict[str, Any]:
    """
    Enrich condition recommendations with detailed reagent information.
    
    Args:
        conditions: Conditions dictionary from recommendation engine
        
    Returns:
        Enriched conditions with reagent details
    """
    enriched = dict(conditions)
    
    # Map condition keys to reagent types
    reagent_type_map = {
        'catalyst': 'metal_catalyst',
        'pd_source': 'metal_catalyst',
        'cu_source': 'metal_catalyst',
        'ni_source': 'metal_catalyst',
        'metal_source': 'metal_catalyst',
        'ligand': 'ligand',
        'ligands': 'ligand',
        'base': 'base',
        'solvent': 'solvent',
        'solvents': 'solvent',
        'oxidant': 'oxidant',
        'reductant': 'reductant',
        'additive': 'additive',
        'acid': 'acid',
    }
    
    # Enrich each condition category
    for key, value in conditions.items():
        reagent_type = reagent_type_map.get(key.lower())
        
        if not reagent_type:
            continue
        
        # Handle list of options
        if isinstance(value, list):
            enriched_list = []
            for item in value:
                if isinstance(item, str):
                    enriched_list.append(enrich_reagent_info(item, reagent_type))
                else:
                    enriched_list.append(item)
            enriched[f"{key}_details"] = enriched_list
        
        # Handle single string value
        elif isinstance(value, str):
            enriched[f"{key}_details"] = enrich_reagent_info(value, reagent_type)
        
        # Handle dict with name key
        elif isinstance(value, dict) and 'name' in value:
            enriched[f"{key}_details"] = enrich_reagent_info(value['name'], reagent_type)
    
    return enriched


def format_reagent_for_display(reagent_info: Dict[str, Any], compact: bool = False) -> str:
    """
    Format reagent information for display.
    
    Args:
        reagent_info: Enriched reagent dictionary
        compact: If True, return compact single-line format
        
    Returns:
        Formatted string
    """
    if not reagent_info.get('found'):
        return reagent_info.get('name', 'Unknown')
    
    parts = [reagent_info['name']]
    
    if reagent_info.get('abbreviation'):
        parts.append(f"({reagent_info['abbreviation']})")
    
    if reagent_info.get('cas'):
        parts.append(f"[CAS: {reagent_info['cas']}]")
    
    if compact:
        return ' '.join(parts)
    
    # Full format with additional details
    lines = [' '.join(parts)]
    
    if reagent_info.get('smiles'):
        lines.append(f"  SMILES: {reagent_info['smiles']}")
    
    if reagent_info.get('aliases') and len(reagent_info['aliases']) > 0:
        lines.append(f"  Also known as: {', '.join(reagent_info['aliases'][:3])}")
    
    return '\n'.join(lines)


def get_all_reagent_types() -> List[str]:
    """Get list of all available reagent database types."""
    reagent_dir = get_data_dir() / "reagent_db"
    if not reagent_dir.exists():
        return []
    
    types = []
    for file in reagent_dir.glob("*.json"):
        if file.stem not in ['not_determined_reagents']:
            types.append(file.stem)
    
    return sorted(types)


def get_all_reagents_by_type(reagent_type: str) -> List[Dict[str, Any]]:
    """
    Get all reagents of a specific type.
    
    Args:
        reagent_type: Type of reagent (e.g., 'ligand', 'base', 'solvent', 'metal_catalyst', 'additive')
        
    Returns:
        List of all reagent dictionaries in that database
        
    Example:
        >>> solvents = get_all_reagents_by_type('solvent')
        >>> print(f"Found {len(solvents)} solvents")
        >>> for s in solvents[:5]:
        ...     print(f"  - {s.get('name')} (CAS: {s.get('cas')})")
    """
    return load_reagent_database(reagent_type)


def count_reagents_by_type(reagent_type: str) -> int:
    """
    Count reagents in a specific database.
    
    Args:
        reagent_type: Type of reagent database
        
    Returns:
        Number of reagents in that database
    """
    db = load_reagent_database(reagent_type)
    return len(db)


def is_reagent_in_database(name: str, reagent_type: str) -> bool:
    """
    Check if a reagent exists in the specified database.
    
    Args:
        name: Reagent name to check
        reagent_type: Type of reagent database to search
        
    Returns:
        True if reagent is found, False otherwise
    """
    if not name:
        return False
    
    reagent = find_reagent(name, reagent_type)
    return reagent is not None


def check_precedent_reagents_in_database(precedent: Dict[str, Any]) -> Dict[str, Any]:
    """
    Check if all reagents in a precedent can be found in the reagent database.
    
    Args:
        precedent: Precedent dictionary with catalytic_system, reagents, and solvents
        
    Returns:
        Dict with:
            - all_found (bool): True if all reagents are in database
            - missing (List[Dict]): List of missing reagents with name and type
            - found_count (int): Number of reagents found
            - total_count (int): Total number of reagents
    """
    missing = []
    found_count = 0
    total_count = 0
    
    # Check catalytic system (metal catalysts and ligands)
    for cat in precedent.get("catalytic_system", []):
        name = cat.get("name")
        if not name:
            continue
            
        total_count += 1
        
        # Try metal_catalyst first, then ligand
        found = (is_reagent_in_database(name, "metal_catalyst") or 
                is_reagent_in_database(name, "ligand"))
        
        if found:
            found_count += 1
        else:
            missing.append({"name": name, "type": "catalyst/ligand"})
    
    # Check reagents (bases, additives, etc.)
    for reagent in precedent.get("reagents", []):
        name = reagent.get("name")
        role = reagent.get("role", "").upper()
        if not name:
            continue
            
        total_count += 1
        
        # Map role to database type
        role_to_type = {
            "BASE": "base",
            "ADDITIVE": "additive",
            "ACID": "acid",
            "OXIDANT": "oxidant",
            "REDUCTANT": "reductant",
        }
        
        reagent_type = role_to_type.get(role, "additive")
        found = is_reagent_in_database(name, reagent_type)
        
        if found:
            found_count += 1
        else:
            missing.append({"name": name, "type": reagent_type})
    
    # Check solvents
    for solvent in precedent.get("solvents", []):
        name = solvent.get("name")
        if not name:
            continue
            
        total_count += 1
        found = is_reagent_in_database(name, "solvent")
        
        if found:
            found_count += 1
        else:
            missing.append({"name": name, "type": "solvent"})
    
    return {
        "all_found": len(missing) == 0 and total_count > 0,
        "missing": missing,
        "found_count": found_count,
        "total_count": total_count,
    }


def filter_precedents_by_database_availability(
    precedents: List[Dict[str, Any]],
    require_all_in_database: bool = True
) -> List[Dict[str, Any]]:
    """
    Filter precedents to only include those where reagents are in the database.
    
    Args:
        precedents: List of precedent dictionaries
        require_all_in_database: If True, require ALL reagents to be in database.
                                If False, keep precedents with at least some reagents found.
        
    Returns:
        Filtered list of precedents
    """
    filtered = []
    
    for prec in precedents:
        check_result = check_precedent_reagents_in_database(prec)
        
        if require_all_in_database:
            # Only keep if all reagents found
            if check_result["all_found"]:
                filtered.append(prec)
        else:
            # Keep if at least some reagents found
            if check_result["found_count"] > 0:
                filtered.append(prec)
    
    return filtered


# Preload common databases on module import (optional)
def preload_common_databases():
    """Preload commonly used reagent databases for faster lookup."""
    common_types = ['ligand', 'base', 'solvent', 'metal_catalyst']
    for reagent_type in common_types:
        try:
            load_reagent_database(reagent_type)
        except Exception:
            pass


if __name__ == "__main__":
    # Test the module
    print("Testing reagent lookup...")
    
    # Test ligand lookup
    print("\n1. Looking up 'PPh3':")
    info = enrich_reagent_info("PPh3", "ligand")
    print(format_reagent_for_display(info))
    
    # Test base lookup
    print("\n2. Looking up 'K2CO3':")
    info = enrich_reagent_info("K2CO3", "base")
    print(format_reagent_for_display(info))
    
    # Test metal catalyst
    print("\n3. Looking up 'Pd(PPh3)4':")
    info = enrich_reagent_info("Pd(PPh3)4", "metal_catalyst")
    print(format_reagent_for_display(info))
    
    # Test condition enrichment
    print("\n4. Enriching full conditions:")
    test_conditions = {
        "catalyst": "Pd(PPh3)4",
        "ligands": ["PPh3", "XPhos"],
        "base": "K2CO3",
        "solvent": ["DMF", "Toluene"]
    }
    enriched = enrich_conditions(test_conditions)
    print(json.dumps(enriched, indent=2, default=str))
