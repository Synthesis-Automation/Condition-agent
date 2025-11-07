"""
Direct test of SNAr condition recommendation without LangChain wrapper.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.detection import detect_reaction
from chemtools.rule.engine import RuleEngine
from chemtools.featurizers.molecular import featurize

print("=" * 70)
print("Direct SNAr Recommendation Test")
print("=" * 70)
print()

# User's SNAr reaction
reaction_smiles = "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1"
reactants = ["Clc1nc(Cl)nc(Cl)n1", "CO"]

print(f"Reaction: {reaction_smiles}")
print()

# Step 1: Detect reaction family
print("Step 1: Detecting reaction family...")
detection = detect_reaction(reaction_smiles)
print(f"  Family: {detection.get('family', 'Unknown')}")
print(f"  Confidence: {detection.get('confidence', 0)}")
print()

# Step 2: Load SNAr rules
print("Step 2: Loading SNAr rule database...")
try:
    from chemtools.rule.database import RuleDatabase
    database = RuleDatabase.from_file("data/rule_db/SNAr_db.json")
    engine = RuleEngine(database)
    print(f"  ✓ Loaded {len(database.base_rules)} base rules")
    print(f"  ✓ Loaded {len(database.modifiers)} modifiers")
except Exception as e:
    print(f"  ✗ Error: {e}")
    sys.exit(1)
print()

# Step 3: Featurize the reaction (for inspection)
print("Step 3: Featurizing reaction (for inspection)...")
features = featurize(reactants[0], reactants[1])
print(f"  Generated {len(features)} features")
print(f"  Sample features: {list(features.keys())[:10]}")
print()

# Step 4: Apply rules (using SMILES, not features dict)
print("Step 4: Applying SNAr rules...")
try:
    recommendation = engine.recommend(reaction_smiles)
    
    if recommendation.base_rule:
        print(f"  ✓ Matched rule: {recommendation.base_rule.name}")
        print(f"    Confidence: {recommendation.base_rule.confidence}")
        print()
        print("  Recommended conditions:")
        for key, value in recommendation.base_rule.conditions.items():
            print(f"    {key}: {value}")
        print()
        
        if recommendation.modifiers:
            print("  Applied modifiers:")
            for modifier in recommendation.modifiers:
                print(f"    - {modifier.suggestion}")
        print()
        
        print("=" * 70)
        print("✅ SUCCESS: SNAr rules applied correctly!")
        print("=" * 70)
    else:
        print("  ✗ No rule matched")
        print()
        print("=" * 70)
        print("⚠️  WARNING: Using default conditions")
        print("=" * 70)
        
except Exception as e:
    print(f"  ✗ Error: {e}")
    import traceback
    traceback.print_exc()
    print()
    print("=" * 70)
    print("❌ FAILED: Could not apply SNAr rules")
    print("=" * 70)
    sys.exit(1)
