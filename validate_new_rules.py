"""
Validate new rule database files
"""
import json
from pathlib import Path
import sys

def validate_rule_file(filepath):
    """Validate a rule database JSON file"""
    print(f"\n{'='*70}")
    print(f"Validating: {filepath.name}")
    print(f"{'='*70}\n")
    
    issues = []
    warnings = []
    
    # 1. Check file exists and is valid JSON
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"❌ ERROR: File not found: {filepath}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ ERROR: Invalid JSON: {e}")
        return False
    
    print("✓ Valid JSON format")
    
    # 2. Check required top-level fields
    required_fields = ["name", "reaction_type", "version", "applies_if"]
    for field in required_fields:
        if field not in data:
            issues.append(f"Missing required field: '{field}'")
        else:
            print(f"✓ Has '{field}': {data[field] if field != 'applies_if' else '...'}")
    
    # 3. Check applies_if structure
    if "applies_if" in data:
        applies = data["applies_if"]
        if not isinstance(applies, dict):
            issues.append("'applies_if' must be a dictionary")
        else:
            if "all" not in applies and "any" not in applies:
                warnings.append("'applies_if' should have 'all' or 'any' keys")
            
            # Check for valid token names
            for key in ["all", "any"]:
                if key in applies:
                    tokens = applies[key]
                    if not isinstance(tokens, list):
                        issues.append(f"'applies_if.{key}' must be a list")
                    else:
                        print(f"  ✓ applies_if.{key}: {len(tokens)} conditions")
    
    # 4. Check default_rule
    if "default_rule" not in data:
        warnings.append("No 'default_rule' defined (optional but recommended)")
    else:
        default = data["default_rule"]
        if "conditions" not in default:
            issues.append("'default_rule' missing 'conditions'")
        else:
            print(f"✓ Has default_rule with {len(default['conditions'])} conditions")
    
    # 5. Check base_rules
    if "base_rules" not in data:
        issues.append("Missing 'base_rules' array")
    else:
        rules = data["base_rules"]
        if not isinstance(rules, list):
            issues.append("'base_rules' must be an array")
        else:
            print(f"✓ Has {len(rules)} base_rules")
            
            # Validate each rule
            for i, rule in enumerate(rules):
                rule_id = rule.get("id", f"rule_{i}")
                
                if "name" not in rule:
                    issues.append(f"base_rules[{i}] ({rule_id}): missing 'name'")
                
                if "conditions" not in rule:
                    issues.append(f"base_rules[{i}] ({rule_id}): missing 'conditions'")
                elif not isinstance(rule["conditions"], dict):
                    issues.append(f"base_rules[{i}] ({rule_id}): 'conditions' must be a dict")
                
                if "reactant_features" not in rule:
                    warnings.append(f"base_rules[{i}] ({rule_id}): no 'reactant_features' (will always match)")
                else:
                    features = rule["reactant_features"]
                    if not isinstance(features, dict):
                        issues.append(f"base_rules[{i}] ({rule_id}): 'reactant_features' must be a dict")
                    else:
                        # Check for valid logic keys
                        valid_keys = {"all", "any", "and", "not"}
                        found_keys = set(features.keys()) & valid_keys
                        if not found_keys:
                            warnings.append(f"base_rules[{i}] ({rule_id}): 'reactant_features' has no logic keys (all/any/and/not)")
    
    # 6. Check modifiers
    if "modifiers" not in data:
        warnings.append("No 'modifiers' array (optional but recommended)")
    else:
        mods = data["modifiers"]
        if not isinstance(mods, list):
            issues.append("'modifiers' must be an array")
        else:
            print(f"✓ Has {len(mods)} modifiers")
            
            # Validate each modifier
            for i, mod in enumerate(mods):
                mod_id = mod.get("id", f"mod_{i}")
                
                if "when" not in mod:
                    issues.append(f"modifiers[{i}] ({mod_id}): missing 'when'")
                elif not isinstance(mod["when"], list):
                    issues.append(f"modifiers[{i}] ({mod_id}): 'when' must be a list")
                
                # Check for suggestion field (could be "suggest" or "suggestion")
                if "suggest" not in mod and "suggestion" not in mod:
                    issues.append(f"modifiers[{i}] ({mod_id}): missing 'suggest' or 'suggestion'")
    
    # 7. Print summary
    print(f"\n{'-'*70}")
    print("VALIDATION SUMMARY:")
    print(f"{'-'*70}")
    
    if not issues and not warnings:
        print("✅ ALL CHECKS PASSED - File is compatible with the system!")
        return True
    
    if warnings:
        print(f"\n⚠️  {len(warnings)} Warning(s):")
        for w in warnings:
            print(f"  • {w}")
    
    if issues:
        print(f"\n❌ {len(issues)} Error(s) - MUST BE FIXED:")
        for e in issues:
            print(f"  • {e}")
        return False
    
    print("\n✅ No critical errors - File should work (with warnings)")
    return True


def main():
    # Check the two new files
    files = [
        Path(r"c:\Git-softwares\Condition-agent\data\rule_db\reductive_amination_db.json"),
        Path(r"c:\Users\xubar\Downloads\SNAr_rule.json")
    ]
    
    results = {}
    for filepath in files:
        if filepath.exists():
            results[filepath.name] = validate_rule_file(filepath)
        else:
            print(f"\n❌ File not found: {filepath}")
            results[filepath.name] = False
    
    # Final summary
    print(f"\n{'='*70}")
    print("FINAL RESULTS:")
    print(f"{'='*70}")
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {name}")
    
    # Exit with error code if any failed
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
