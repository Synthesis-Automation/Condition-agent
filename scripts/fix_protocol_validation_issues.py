"""
Fix validation issues in migrated v2.0 protocol files.

Common issues fixed:
1. Invalid chemical roles (nucleophile, electrophile → starting_material)
2. Other schema compliance issues

Usage:
    python scripts/fix_protocol_validation_issues.py
"""

import json
from pathlib import Path
from typing import Any, Dict, List


# Map non-standard roles to standard ones
ROLE_MAPPING = {
    "nucleophile": "starting_material",
    "electrophile": "starting_material",
    "boronic acid": "starting_material",
    "coupling partner": "starting_material",
    "substrate": "starting_material",
    "reductant": "reagent",
    "oxidant": "reagent",
    "acid": "additive",
    "photocatalyst": "catalyst",
    "activator": "reagent",
    "photoredox catalyst": "catalyst",
}


def fix_chemical_roles(data: Dict[str, Any]) -> bool:
    """Fix non-standard chemical roles in reaction_setup."""
    modified = False
    
    for setup in data.get("reaction_setup", []):
        for chemical in setup.get("chemicals", []):
            role = chemical.get("role", "")
            
            # Check if role needs mapping
            if role in ROLE_MAPPING:
                chemical["role"] = ROLE_MAPPING[role]
                modified = True
                print(f"    Fixed role: {role} → {chemical['role']}")
    
    return modified


def fix_null_conditions(data: Dict[str, Any]) -> bool:
    """Remove or fix condition entries with null temperature or atmosphere."""
    modified = False
    
    for setup in data.get("reaction_setup", []):
        conditions = setup.get("conditions", [])
        
        # Fix conditions with null values
        for cond in conditions:
            # Remove null atmosphere (optional field)
            if "atmosphere" in cond and cond["atmosphere"] is None:
                del cond["atmosphere"]
                modified = True
                print(f"    Removed null atmosphere from condition")
        
        # Filter out conditions with null temperature (required field)
        original_count = len(conditions)
        setup["conditions"] = [
            cond for cond in conditions
            if cond.get("temperature_C") is not None
        ]
        
        if len(setup["conditions"]) < original_count:
            removed = original_count - len(setup["conditions"])
            print(f"    Removed {removed} condition(s) with null temperature")
            modified = True
    
    return modified


def fix_workup_strings(data: Dict[str, Any]) -> bool:
    """Convert string workup steps to structured objects."""
    modified = False
    
    for workup_section in data.get("workup_and_purification", []):
        # Fix quench steps
        if "quench" in workup_section:
            quench_steps = workup_section["quench"]
            fixed_steps = []
            had_strings = False
            
            for step in quench_steps:
                if isinstance(step, str):
                    # Convert string to object
                    fixed_steps.append({
                        "reagent": "Quench reagent",
                        "details": step
                    })
                    had_strings = True
                elif isinstance(step, dict):
                    fixed_steps.append(step)
            
            if had_strings:
                workup_section["quench"] = fixed_steps
                print(f"    Fixed string quench steps → structured objects")
                modified = True
        
        # Fix workup steps
        if "workup" in workup_section:
            workup_steps = workup_section["workup"]
            fixed_steps = []
            had_strings = False
            
            for step in workup_steps:
                if isinstance(step, str):
                    # Convert string to object
                    fixed_steps.append({
                        "step": "Processing",
                        "details": step
                    })
                    had_strings = True
                elif isinstance(step, dict):
                    fixed_steps.append(step)
            
            if had_strings:
                workup_section["workup"] = fixed_steps
                print(f"    Fixed string workup steps → structured objects")
                modified = True
        
        # Fix purification steps
        if "purification" in workup_section:
            purif_steps = workup_section["purification"]
            fixed_steps = []
            had_strings = False
            
            for step in purif_steps:
                if isinstance(step, str):
                    # Convert string to object
                    fixed_steps.append({
                        "method": "Purification",
                        "details": step
                    })
                    had_strings = True
                elif isinstance(step, dict):
                    fixed_steps.append(step)
            
            if had_strings:
                workup_section["purification"] = fixed_steps
                print(f"    Fixed string purification steps → structured objects")
                modified = True
    
    return modified


def fix_protocol_file(file_path: Path) -> bool:
    """Fix validation issues in a single protocol file."""
    
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    # Handle array format
    if isinstance(data, list):
        any_modified = False
        for idx, protocol in enumerate(data):
            print(f"  Protocol {idx}:")
            modified1 = fix_chemical_roles(protocol)
            modified2 = fix_null_conditions(protocol)
            modified3 = fix_workup_strings(protocol)
            if modified1 or modified2 or modified3:
                any_modified = True
        
        if any_modified:
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            return True
    else:
        # Single protocol
        modified1 = fix_chemical_roles(data)
        modified2 = fix_null_conditions(data)
        modified3 = fix_workup_strings(data)
        if modified1 or modified2 or modified3:
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            return True
    
    return False


def fix_all_protocols():
    """Fix validation issues in all protocol files."""
    
    protocol_dir = Path("data/protocol_db_v2")
    
    # Find all JSON files
    protocol_files = [
        f for f in protocol_dir.glob("*.json")
        if not f.name.startswith(".")
    ]
    
    print(f"Checking {len(protocol_files)} protocol files for validation issues...\n")
    
    fixed_count = 0
    unchanged_count = 0
    
    for protocol_file in sorted(protocol_files):
        print(f"Checking: {protocol_file.name}")
        
        try:
            modified = fix_protocol_file(protocol_file)
            
            if modified:
                fixed_count += 1
                print(f"  ✅ Fixed and saved\n")
            else:
                unchanged_count += 1
                print(f"  ⏭️  No issues found\n")
        
        except Exception as e:
            print(f"  ❌ Error: {e}\n")
    
    print(f"{'='*60}")
    print(f"Fix Summary:")
    print(f"  ✅ Fixed: {fixed_count} files")
    print(f"  ⏭️  Unchanged: {unchanged_count} files")
    print(f"  📁 Total: {len(protocol_files)} files")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    fix_all_protocols()
