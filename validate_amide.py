#!/usr/bin/env python
"""Quick validation of amide_formation.json structure."""

import json
from pathlib import Path

def validate_amide_formation():
    """Validate amide_formation.json structure."""
    
    print("="*70)
    print("AMIDE_FORMATION.JSON STRUCTURE VALIDATION")
    print("="*70)
    
    # Load file
    db_path = Path("data/rule_db/amide_formation.json")
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
    metadata_fields = ["name", "reaction_type", "version", "evaluation"]
    for field in metadata_fields:
        if field in db:
            value = db[field]
            if len(str(value)) > 60:
                print(f"   ✓ {field}: {str(value)[:60]}...")
            else:
                print(f"   ✓ {field}: {value}")
    
    # 2. Check applies_if structure
    print("\n2. applies_if validation:")
    if "applies_if" in db:
        applies_if = db["applies_if"]
        if "all" in applies_if:
            features = applies_if["all"]
            print(f"   ✓ Requires ALL of: {features}")
        if "any" in applies_if:
            features = applies_if["any"]
            print(f"   ✓ Requires ANY of: {features}")
        if not ("all" in applies_if or "any" in applies_if):
            warnings.append("applies_if has no 'all' or 'any' key")
    
    # 3. Check default_rule
    print("\n3. default_rule validation:")
    if "default_rule" in db:
        default_rule = db["default_rule"]
        if "conditions" in default_rule:
            print(f"   ✓ Has conditions with {len(default_rule['conditions'])} parameters")
            for key in list(default_rule["conditions"].keys())[:5]:
                value = default_rule["conditions"][key]
                val_str = str(value)[:60] + "..." if len(str(value)) > 60 else str(value)
                print(f"      • {key}: {val_str}")
        else:
            issues.append("default_rule missing 'conditions'")
        
        if "id" in default_rule:
            print(f"   ✓ ID: {default_rule['id']}")
        if "description" in default_rule:
            print(f"   ✓ Description: {default_rule['description'][:60]}...")
    
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
                print(f"   ✓ Matches when ANY: {features['any'][:3]}...")
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
        if i > 5:  # Only show first 5
            print(f"\n   ... and {len(db.get('modifiers', [])) - 5} more modifiers")
            break
        print(f"\n   Modifier {i}:")
        
        # Check id
        if "id" in mod:
            print(f"   ✓ id: {mod['id']}")
        
        # Check when
        if "when" in mod:
            conditions = mod["when"]
            print(f"   ✓ when: {conditions}")
        else:
            issues.append(f"modifiers[{i-1}] missing 'when'")
        
        # Check suggest
        if "suggest" in mod:
            print(f"   ✓ suggest: {mod['suggest'][:50]}...")
        else:
            issues.append(f"modifiers[{i-1}] missing 'suggest'")
    
    # 6. Check feature_detection section
    print("\n6. feature_detection section:")
    if "feature_detection" in db:
        print(f"   ✓ feature_detection present")
        fd = db["feature_detection"]
        if "detectors" in fd:
            print(f"   ✓ {len(fd['detectors'])} detectors defined")
            # Check first few detectors
            for i, detector in enumerate(fd["detectors"][:3], 1):
                if "id" in detector:
                    print(f"      • {detector['id']}")
        if "notes" in fd:
            print(f"   ✓ notes: {fd['notes'][:50]}...")
    else:
        warnings.append("No feature_detection section (optional but recommended)")
    
    # 7. Check applies_if_smarts section
    print("\n7. applies_if_smarts section:")
    if "applies_if_smarts" in db:
        print(f"   ✓ applies_if_smarts present")
        ais = db["applies_if_smarts"]
        if "all" in ais:
            print(f"   ✓ {len(ais['all'])} required features")
        if "any" in ais:
            print(f"   ✓ {len(ais['any'])} alternative features")
    else:
        warnings.append("No applies_if_smarts section (optional)")
    
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
        print(f"\n⚠️  WARNINGS ({len(warnings)}):")
        for warning in warnings:
            print(f"   - {warning}")
    else:
        print("\n✅ No warnings")
    
    if not issues:
        print("\n🎉 Database structure is valid!")
    
    return len(issues) == 0

if __name__ == "__main__":
    success = validate_amide_formation()
    exit(0 if success else 1)
