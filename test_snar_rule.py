"""
Test SNAr rule database integration
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.rule.engine import RuleEngine

# Test 1: Can we load the SNAr rule file directly?
print("="*70)
print("Test 1: Loading SNAr rule database directly")
print("="*70)
try:
    engine = RuleEngine.from_file("data/rule_db/SNAr_db.json")
    print("✓ Successfully loaded SNAr_db.json")
    print(f"  Name: {engine.database.metadata.get('name')}")
    print(f"  Type: {engine.database.metadata.get('reaction_type')}")
    print(f"  Base rules: {len(engine.database.base_rules)}")
    print(f"  Modifiers: {len(engine.database.modifiers)}")
except Exception as e:
    print(f"✗ Failed to load: {e}")
    sys.exit(1)

# Test 2: Check the wrapper mapping
print("\n" + "="*70)
print("Test 2: Testing family-to-database mapping")
print("="*70)

from chem_assistant.chemtools_wrapper import _FAMILY_TO_RULE_DB

test_families = ["snar", "s_nar", "SNAr", "S_NAr", "aromatic_nucleophilic_substitution"]
for family in test_families:
    normalized = family.lower().replace("-", "_").replace(" ", "_")
    mapped_db = _FAMILY_TO_RULE_DB.get(normalized, "NOT FOUND")
    print(f"  {family:40s} -> {mapped_db}")

# Test 3: Test feature detection for SNAr
print("\n" + "="*70)
print("Test 3: Testing feature detection for SNAr reaction")
print("="*70)

reaction_smiles = "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1"
print(f"Reaction: {reaction_smiles}")

from chemtools.featurizers.molecular import featurize

try:
    # Parse the reaction
    parts = reaction_smiles.split(">>")
    reactants = parts[0].split(".")
    
    print(f"\nElectrophile: {reactants[0]}")
    print(f"Nucleophile: {reactants[1]}")
    
    # Featurize the reaction
    features = featurize(reactants[0], reactants[1])
    
    snar_tokens = [k for k in features.keys() if 'snar' in k.lower() or 'aromatic' in k.lower()]
    print(f"\nSNAr-related features found: {len(snar_tokens)}")
    for token in snar_tokens[:10]:  # Show first 10
        print(f"  {token}: {features[token]}")
    
    nuc_tokens = [k for k in features.keys() if features[k] and k.endswith('_present')]
    print(f"\nNucleophile features found: {len(nuc_tokens)}")
    for token in nuc_tokens[:10]:  # Show first 10
        print(f"  {token}: {features[token]}")
    
except Exception as e:
    print(f"✗ Feature extraction failed: {e}")
    import traceback
    traceback.print_exc()

# Test 4: Try to get recommendation using the engine
print("\n" + "="*70)
print("Test 4: Getting recommendation from SNAr engine")
print("="*70)

try:
    # Combine features for the reaction (already computed above)
    all_features = features
    
    # Check if SNAr applies
    applies = engine.database.check_applies(all_features)
    print(f"SNAr database applies: {applies}")
    
    if applies:
        # Find matching rule
        matched_rule = engine.database.find_matching_rule(all_features)
        if matched_rule:
            print(f"\n✓ Matched rule: {matched_rule.name}")
            print(f"  Conditions:")
            for key, value in matched_rule.conditions.items():
                print(f"    {key}: {value}")
        else:
            print("  No specific rule matched, would use default")
    
except Exception as e:
    print(f"✗ Recommendation failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("All tests completed!")
print("="*70)
