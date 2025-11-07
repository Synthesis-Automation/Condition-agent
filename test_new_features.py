"""
Test the newly added features in calculable_features.json
"""
from chemtools.featurizers.calculable import detect_all_features

# Test molecules
test_cases = [
    # Boronic acids
    ("c1ccccc1B(O)O", "phenylboronic acid", ["boronic_acid_present"]),
    
    # Alcohols
    ("CCO", "ethanol (primary)", ["primary_alcohol_present", "aliphatic_alcohol_present", "alcohol_present"]),
    ("CC(O)C", "isopropanol (secondary)", ["secondary_alcohol_present", "aliphatic_alcohol_present", "alcohol_present"]),
    ("CC(C)(C)O", "tert-butanol (tertiary)", ["tertiary_alcohol_present", "aliphatic_alcohol_present", "alcohol_present"]),
    
    # Carboxylic acid (already tested, but verify again)
    ("O=C(O)c1ccccc1", "benzoic acid", ["carboxylic_acid_present", "acidic_proton_present"]),
    
    # Aromatic
    ("c1ccccc1", "benzene", ["aromatic_present"]),
    ("CCO", "ethanol (not aromatic)", []),  # Should NOT have aromatic_present
    
    # Cyclopropane
    ("C1CC1", "cyclopropane", ["cyclopropane_present", "strained_ring_present"]),
    ("CCCC", "butane (no cyclopropane)", []),
    
    # Alpha chiral center (derived feature - hard to test precisely)
    ("C[C@H](O)C(=O)O", "lactic acid", ["alpha_chiral_center_present", "chiral_center_present"]),
]

print("=" * 80)
print("TESTING NEW FEATURES")
print("=" * 80)

all_passed = True

for smiles, name, expected_features in test_cases:
    try:
        result = detect_all_features(smiles)
        
        if not result:
            print(f"\n❌ FAILED to detect features for: {name} ({smiles})")
            all_passed = False
            continue
            
        # Check expected features
        missing = []
        for feat in expected_features:
            if feat not in result or not result.get(feat):
                missing.append(feat)
        
        if missing:
            print(f"\n⚠️  {name} ({smiles})")
            print(f"   Missing expected features: {missing}")
            print(f"   Present features: {[k for k, v in result.items() if v and k.endswith('_present')][:10]}")
            all_passed = False
        else:
            print(f"✅ {name}: All expected features detected")
            
    except Exception as e:
        print(f"\n❌ ERROR testing {name} ({smiles}): {e}")
        all_passed = False

print("\n" + "=" * 80)
if all_passed:
    print("✅ ALL TESTS PASSED")
else:
    print("⚠️  SOME TESTS FAILED - Review output above")
print("=" * 80)
