"""
Final validation: Complete integration test for SNAr and Reductive Amination rule databases.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.detection import detect_reaction
from chemtools.rule.engine import RuleEngine
from chemtools.rule.database import RuleDatabase

print("=" * 80)
print("FINAL VALIDATION: SNAr and Reductive Amination Rule Database Integration")
print("=" * 80)
print()

# Test cases
test_cases = [
    {
        "name": "SNAr - Heteroaryl Chloride + Alcohol",
        "reaction": "Clc1nc(Cl)nc(Cl)n1.CO>>COc1nc(Cl)nc(Cl)n1",
        "expected_family": "snar",
        "database_file": "data/rule_db/SNAr_db.json"
    },
    {
        "name": "Reductive Amination - Benzaldehyde + Aniline",
        "reaction": "O=Cc1ccccc1.Nc1ccccc1>>c1ccccc1CNc1ccccc1",
        "expected_family": "reductive_amination",
        "database_file": "data/rule_db/reductive_amination_db.json"
    }
]

all_passed = True

for i, test in enumerate(test_cases, 1):
    print(f"Test {i}: {test['name']}")
    print("-" * 80)
    print(f"Reaction: {test['reaction']}")
    print()
    
    # Step 1: Detection
    print("  [1] Detecting reaction family...")
    try:
        detection = detect_reaction(test['reaction'])
        detected_family = detection.get('family', 'Unknown')
        confidence = detection.get('confidence', 0)
        
        if detected_family == test['expected_family']:
            print(f"      ✓ Detected as '{detected_family}' (confidence: {confidence})")
        else:
            print(f"      ✗ Expected '{test['expected_family']}', got '{detected_family}'")
            all_passed = False
            print()
            continue
    except Exception as e:
        print(f"      ✗ Detection failed: {e}")
        all_passed = False
        print()
        continue
    
    # Step 2: Load rule database
    print(f"  [2] Loading rule database: {test['database_file']}")
    try:
        database = RuleDatabase.from_file(test['database_file'])
        engine = RuleEngine(database)
        print(f"      ✓ Loaded {len(database.base_rules)} base rules, {len(database.modifiers)} modifiers")
    except Exception as e:
        print(f"      ✗ Failed to load database: {e}")
        all_passed = False
        print()
        continue
    
    # Step 3: Apply rules
    print("  [3] Applying rules...")
    try:
        recommendation = engine.recommend(test['reaction'])
        
        if recommendation.base_rule:
            print(f"      ✓ Matched rule: '{recommendation.base_rule.name}'")
            print(f"        Confidence: {recommendation.base_rule.confidence}")
            print(f"        Conditions:")
            for key, value in list(recommendation.base_rule.conditions.items())[:3]:
                # Show first 3 conditions
                value_str = str(value)[:60] + "..." if len(str(value)) > 60 else str(value)
                print(f"          • {key}: {value_str}")
        else:
            print(f"      ✗ No matching rule found")
            all_passed = False
    except Exception as e:
        print(f"      ✗ Rule application failed: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False
    
    print()

print("=" * 80)
if all_passed:
    print("✅ ALL TESTS PASSED - Integration complete and functional!")
else:
    print("❌ SOME TESTS FAILED - Review errors above")
print("=" * 80)

sys.exit(0 if all_passed else 1)
