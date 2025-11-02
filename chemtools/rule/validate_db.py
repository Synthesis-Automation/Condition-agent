#!/usr/bin/env python
"""
Rule Database Structure Validator
==================================

Comprehensive validation of rule database JSON files.
"""

import json
from pathlib import Path
from typing import Dict, List, Any

def validate_suzuki_json():
    """Validate suzuki.json for logical consistency."""
    
    print("="*70)
    print("SUZUKI.JSON STRUCTURE VALIDATION")
    print("="*70)
    
    # Load file
    db_path = Path("data/rule_db/suzuki.json")
    with open(db_path, encoding='utf-8') as f:
        db = json.load(f)
    
    issues = []
    warnings = []
    
    # 1. Check required top-level fields
    print("\n1. Top-level structure:")
    required_fields = ["applies_if", "default_rule", "base_rules", "modifiers"]
    for field in required_fields:
        if field in db:
            print(f"   ✓ {field}")
        else:
            issues.append(f"Missing required field: {field}")
            print(f"   ✗ {field} - MISSING")
    
    # Optional metadata
    metadata_fields = ["name", "reaction_type", "version"]
    for field in metadata_fields:
        if field in db:
            print(f"   ✓ {field}: {db[field]}")
    
    # 2. Check applies_if structure
    print("\n2. applies_if validation:")
    if "applies_if" in db:
        applies_if = db["applies_if"]
        if "all" in applies_if:
            features = applies_if["all"]
            print(f"   ✓ Requires ALL of: {features}")
            if not isinstance(features, list):
                issues.append("applies_if['all'] must be a list")
            elif len(features) == 0:
                warnings.append("applies_if['all'] is empty")
        elif "any" in applies_if:
            features = applies_if["any"]
            print(f"   ✓ Requires ANY of: {features}")
        else:
            warnings.append("applies_if has no 'all' or 'any' key")
    
    # 3. Check default_rule
    print("\n3. default_rule validation:")
    if "default_rule" in db:
        default_rule = db["default_rule"]
        if "conditions" in default_rule:
            print(f"   ✓ Has conditions with {len(default_rule['conditions'])} parameters")
            for key, value in default_rule["conditions"].items():
                print(f"      • {key}: {value}")
        else:
            issues.append("default_rule missing 'conditions'")
        
        if "id" in default_rule:
            print(f"   ✓ ID: {default_rule['id']}")
        if "description" in default_rule:
            print(f"   ✓ Description: {default_rule['description'][:50]}...")
    
    # 4. Check base_rules
    print(f"\n4. base_rules validation ({len(db.get('base_rules', []))} rules):")
    for i, rule in enumerate(db.get("base_rules", []), 1):
        print(f"\n   Rule {i}:")
        
        # Check name
        if "name" in rule:
            print(f"   ✓ name: {rule['name']}")
        else:
            warnings.append(f"base_rules[{i-1}] missing 'name' field")
        
        # Check id
        if "id" in rule:
            print(f"   ✓ id: {rule['id']}")
        else:
            warnings.append(f"base_rules[{i-1}] missing 'id' field")
        
        # Check reactant_features
        if "reactant_features" in rule:
            features = rule["reactant_features"]
            if "and" in features:
                print(f"   ✓ Matches when ALL: {features['and']}")
            elif "any" in features:
                print(f"   ✓ Matches when ANY: {features['any']}")
            elif "all" in features:
                print(f"   ✓ Matches when ALL: {features['all']}")
            else:
                warnings.append(f"base_rules[{i-1}] reactant_features has no and/any/all")
        else:
            warnings.append(f"base_rules[{i-1}] missing 'reactant_features'")
        
        # Check conditions
        if "conditions" in rule:
            print(f"   ✓ conditions: {len(rule['conditions'])} parameters")
        else:
            issues.append(f"base_rules[{i-1}] missing 'conditions'")
    
    # 5. Check modifiers
    print(f"\n5. modifiers validation ({len(db.get('modifiers', []))} modifiers):")
    for i, mod in enumerate(db.get("modifiers", []), 1):
        print(f"\n   Modifier {i}:")
        
        # Check when
        if "when" in mod:
            conditions = mod["when"]
            print(f"   ✓ when: {conditions}")
            
            # Check for symptom vs feature
            for cond in conditions:
                if cond.startswith("symptom:"):
                    print(f"      → Symptom-based: {cond}")
                else:
                    print(f"      → Feature-based: {cond}")
        else:
            issues.append(f"modifiers[{i-1}] missing 'when'")
        
        # Check suggest
        if "suggest" in mod:
            print(f"   ✓ suggest: {mod['suggest'][:50]}...")
        elif "suggestion" in mod:
            print(f"   ✓ suggestion: {mod['suggestion'][:50]}...")
        else:
            issues.append(f"modifiers[{i-1}] missing 'suggest' or 'suggestion'")
        
        # Check optional fields
        if "rationale" in mod:
            print(f"   ✓ rationale: {mod['rationale'][:50]}...")
        if "id" in mod:
            print(f"   ✓ id: {mod['id']}")
    
    # 6. Check for feature references
    print("\n6. Feature reference check:")
    all_features = set()
    
    # From applies_if
    if "applies_if" in db:
        for key in ["all", "any"]:
            if key in db["applies_if"]:
                all_features.update(db["applies_if"][key])
    
    # From base_rules
    for rule in db.get("base_rules", []):
        if "reactant_features" in rule:
            for key in ["and", "any", "all"]:
                if key in rule["reactant_features"]:
                    all_features.update(rule["reactant_features"][key])
    
    # From modifiers (non-symptom)
    for mod in db.get("modifiers", []):
        if "when" in mod:
            for cond in mod["when"]:
                if not cond.startswith("symptom:"):
                    all_features.add(cond)
    
    print(f"   Total unique features referenced: {len(all_features)}")
    print(f"   Features: {sorted(list(all_features)[:5])}...")
    
    # 7. Logical consistency checks
    print("\n7. Logical consistency:")
    
    # Check if base_rules requirements are subsets of applies_if
    if "applies_if" in db and "all" in db["applies_if"]:
        required_global = set(db["applies_if"]["all"])
        for i, rule in enumerate(db.get("base_rules", []), 1):
            if "reactant_features" in rule and "and" in rule["reactant_features"]:
                rule_features = set(rule["reactant_features"]["and"])
                if not rule_features.issubset(required_global):
                    extra = rule_features - required_global
                    print(f"   ⚠ Rule {i} requires features beyond applies_if: {extra}")
                else:
                    print(f"   ✓ Rule {i} features are consistent with applies_if")
    
    # 8. Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    if issues:
        print(f"\n❌ CRITICAL ISSUES ({len(issues)}):")
        for issue in issues:
            print(f"   - {issue}")
    else:
        print("\n✅ No critical issues found")
    
    if warnings:
        print(f"\n⚠ WARNINGS ({len(warnings)}):")
        for warning in warnings:
            print(f"   - {warning}")
    else:
        print("\n✅ No warnings")
    
    if not issues and not warnings:
        print("\n🎉 Database is structurally sound and logically consistent!")
    
    return len(issues) == 0


if __name__ == "__main__":
    success = validate_suzuki_json()
    exit(0 if success else 1)
