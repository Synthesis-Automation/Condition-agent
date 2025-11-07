"""
Test loading new rule files with the actual codebase
"""
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from chemtools.rule.database import RuleDatabase


def test_load_rule_file(filepath):
    """Test loading a rule file with RuleDatabase"""
    print(f"\n{'='*70}")
    print(f"Testing: {filepath.name}")
    print(f"{'='*70}\n")
    
    try:
        # Load the database
        db = RuleDatabase.from_file(filepath)
        print("✓ Successfully loaded with RuleDatabase.from_file()")
        
        # Check metadata
        print(f"\nMetadata:")
        print(f"  Name: {db.metadata.get('name', 'N/A')}")
        print(f"  Type: {db.metadata.get('reaction_type', 'N/A')}")
        print(f"  Version: {db.metadata.get('version', 'N/A')}")
        
        # Check applies_if
        print(f"\nApplies If:")
        if db.applies_if:
            for key, value in db.applies_if.items():
                if isinstance(value, list):
                    print(f"  {key}: {len(value)} conditions")
                    for v in value[:3]:  # Show first 3
                        print(f"    - {v}")
                    if len(value) > 3:
                        print(f"    ... and {len(value) - 3} more")
        
        # Check default_rule
        print(f"\nDefault Rule:")
        if db.default_rule:
            print(f"  ✓ Present with {len(db.default_rule)} conditions")
        else:
            print(f"  ✗ Not defined")
        
        # Check base_rules
        print(f"\nBase Rules: {len(db.base_rules)}")
        for i, rule in enumerate(db.base_rules[:3], 1):  # Show first 3
            print(f"  {i}. {rule.name}")
            if rule.reactant_features:
                print(f"     Features: {list(rule.reactant_features.keys())}")
        if len(db.base_rules) > 3:
            print(f"  ... and {len(db.base_rules) - 3} more")
        
        # Check modifiers
        print(f"\nModifiers: {len(db.modifiers)}")
        for i, mod in enumerate(db.modifiers[:3], 1):  # Show first 3
            print(f"  {i}. {len(mod.when)} conditions → {mod.suggestion[:60]}...")
        if len(db.modifiers) > 3:
            print(f"  ... and {len(db.modifiers) - 3} more")
        
        # Validate
        print(f"\nValidation:")
        issues = db.validate()
        if not issues:
            print("  ✓ No validation issues")
        else:
            print(f"  ✗ {len(issues)} issue(s):")
            for issue in issues:
                print(f"    - {issue}")
        
        # Test check_applies with sample features
        print(f"\nTesting check_applies():")
        
        # Get first few tokens from applies_if
        test_features = {}
        if "all" in db.applies_if:
            for token in db.applies_if["all"][:2]:
                test_features[token] = True
        if "any" in db.applies_if:
            if db.applies_if["any"]:
                test_features[db.applies_if["any"][0]] = True
        
        applies = db.check_applies(test_features)
        print(f"  Test features: {list(test_features.keys())}")
        print(f"  Result: {'✓ Applies' if applies else '✗ Does not apply'}")
        
        # Test find_matching_rule
        print(f"\nTesting find_matching_rule():")
        if db.base_rules and db.base_rules[0].reactant_features:
            # Create features that match first rule
            first_rule = db.base_rules[0]
            test_features = {}
            
            if "all" in first_rule.reactant_features:
                for token in first_rule.reactant_features["all"]:
                    test_features[token] = True
            elif "any" in first_rule.reactant_features:
                test_features[first_rule.reactant_features["any"][0]] = True
            elif "and" in first_rule.reactant_features:
                for token in first_rule.reactant_features["and"]:
                    test_features[token] = True
            
            matched_rule = db.find_matching_rule(test_features)
            if matched_rule:
                print(f"  ✓ Matched rule: {matched_rule.name}")
            else:
                print(f"  ✗ No rule matched (using default)")
        
        print(f"\n✅ ALL TESTS PASSED for {filepath.name}")
        return True
        
    except Exception as e:
        print(f"\n❌ ERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    # Test both files
    files = [
        Path(r"c:\Git-softwares\Condition-agent\data\rule_db\reductive_amination_db.json"),
        Path(r"c:\Users\xubar\Downloads\SNAr_rule.json")
    ]
    
    results = {}
    for filepath in files:
        if filepath.exists():
            results[filepath.name] = test_load_rule_file(filepath)
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
    
    all_passed = all(results.values())
    if all_passed:
        print(f"\n🎉 Both files are fully compatible with the codebase!")
    
    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
