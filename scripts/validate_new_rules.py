#!/usr/bin/env python
"""
Validate New Rule Database Files
=================================

Script to validate the structure and content of new rule database JSON files.
"""

import json
import sys
from pathlib import Path
from typing import List, Dict, Any

def validate_json_structure(filepath: Path) -> tuple[bool, List[str]]:
    """Validate JSON structure and return (is_valid, errors)."""
    errors = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, [f"JSON parse error: {e}"]
    except FileNotFoundError:
        return False, [f"File not found: {filepath}"]
    except Exception as e:
        return False, [f"Error reading file: {e}"]
    
    # Required top-level fields
    required_fields = ["name", "reaction_type", "version", "applies_if", "default_rule", "base_rules"]
    for field in required_fields:
        if field not in data:
            errors.append(f"Missing required field: '{field}'")
    
    # Validate applies_if structure
    if "applies_if" in data:
        applies_if = data["applies_if"]
        if not isinstance(applies_if, dict):
            errors.append("'applies_if' must be a dictionary")
        else:
            if "all" not in applies_if and "any" not in applies_if:
                errors.append("'applies_if' must contain at least 'all' or 'any' key")
            if "all" in applies_if and not isinstance(applies_if["all"], list):
                errors.append("'applies_if.all' must be a list")
            if "any" in applies_if and not isinstance(applies_if["any"], list):
                errors.append("'applies_if.any' must be a list")
    
    # Validate default_rule
    if "default_rule" in data:
        default_rule = data["default_rule"]
        if not isinstance(default_rule, dict):
            errors.append("'default_rule' must be a dictionary")
        else:
            if "id" not in default_rule:
                errors.append("'default_rule' missing 'id' field")
            if "description" not in default_rule:
                errors.append("'default_rule' missing 'description' field")
            if "conditions" not in default_rule:
                errors.append("'default_rule' missing 'conditions' field")
    
    # Validate base_rules
    if "base_rules" in data:
        base_rules = data["base_rules"]
        if not isinstance(base_rules, list):
            errors.append("'base_rules' must be a list")
        else:
            for idx, rule in enumerate(base_rules):
                if not isinstance(rule, dict):
                    errors.append(f"base_rules[{idx}] must be a dictionary")
                    continue
                
                # Check required fields in each base_rule
                for field in ["name", "id", "description", "conditions"]:
                    if field not in rule:
                        errors.append(f"base_rules[{idx}] missing '{field}' field")
                
                # reactant_features is optional but should be dict if present
                if "reactant_features" in rule and not isinstance(rule["reactant_features"], dict):
                    errors.append(f"base_rules[{idx}].reactant_features must be a dictionary")
    
    # Validate modifiers if present
    if "modifiers" in data:
        modifiers = data["modifiers"]
        if not isinstance(modifiers, list):
            errors.append("'modifiers' must be a list")
        else:
            for idx, mod in enumerate(modifiers):
                if not isinstance(mod, dict):
                    errors.append(f"modifiers[{idx}] must be a dictionary")
                    continue
                
                # Modifiers require 'id', 'when' (list), and 'suggest' (string)
                if "id" not in mod:
                    errors.append(f"modifiers[{idx}] missing 'id' field")
                if "when" not in mod:
                    errors.append(f"modifiers[{idx}] missing 'when' field")
                elif not isinstance(mod["when"], list):
                    errors.append(f"modifiers[{idx}].when must be a list")
                if "suggest" not in mod:
                    errors.append(f"modifiers[{idx}] missing 'suggest' field")
    
    return len(errors) == 0, errors


def check_reaction_type_mapping(filepath: Path, data: Dict[str, Any]) -> List[str]:
    """Check if reaction_type needs to be added to router or detection systems."""
    warnings = []
    
    reaction_type = data.get("reaction_type", "").lower()
    filename = filepath.stem
    
    # Common reaction type variations
    if "c_o" in filename.lower() or "c-o" in reaction_type:
        warnings.append("C-O coupling detected. Ensure 'C_O_coupling' or similar is mapped in router.py")
    elif "rcm" in filename.lower() or "metathesis" in reaction_type:
        warnings.append("RCM detected. Ensure 'ring_closing_metathesis' or 'RCM' is mapped in router.py")
    elif "sonogashira" in filename.lower() or "sonogashira" in reaction_type:
        warnings.append("Sonogashira detected. Ensure 'Sonogashira' family is mapped in router.py")
    
    return warnings


def main():
    """Main validation routine."""
    # Find the three new rule files
    rule_db_dir = Path("data/rule_db")
    
    if not rule_db_dir.exists():
        print(f"Error: Rule database directory not found: {rule_db_dir}")
        sys.exit(1)
    
    files_to_check = [
        rule_db_dir / "C_O_coupling_db.json",
        rule_db_dir / "RCM_db.json",
        rule_db_dir / "sonogashira_db.json"
    ]
    
    print("=" * 70)
    print("Validating New Rule Database Files")
    print("=" * 70)
    
    all_valid = True
    
    for filepath in files_to_check:
        print(f"\n{'='*70}")
        print(f"Checking: {filepath.name}")
        print(f"{'='*70}")
        
        if not filepath.exists():
            print(f"❌ File not found: {filepath}")
            all_valid = False
            continue
        
        # Validate JSON structure
        is_valid, errors = validate_json_structure(filepath)
        
        if is_valid:
            print("✅ JSON structure is valid")
            
            # Load data for additional checks
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            print(f"   Name: {data.get('name', 'N/A')}")
            print(f"   Reaction Type: {data.get('reaction_type', 'N/A')}")
            print(f"   Version: {data.get('version', 'N/A')}")
            print(f"   Base Rules: {len(data.get('base_rules', []))}")
            print(f"   Modifiers: {len(data.get('modifiers', []))}")
            
            # Check for integration warnings
            warnings = check_reaction_type_mapping(filepath, data)
            if warnings:
                print("\n⚠️  Integration Warnings:")
                for warning in warnings:
                    print(f"   - {warning}")
        else:
            print("❌ Validation failed with errors:")
            for error in errors:
                print(f"   - {error}")
            all_valid = False
    
    print(f"\n{'='*70}")
    if all_valid:
        print("✅ All rule files passed validation!")
        print("\nNext steps:")
        print("1. Ensure reaction types are mapped in chemtools/router.py")
        print("2. Update detection logic if needed in chemtools/detection.py")
        print("3. Test with sample reactions using chemtools.rule.cli")
        sys.exit(0)
    else:
        print("❌ Some files failed validation. Please fix the errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
