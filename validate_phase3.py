"""Validate Phase 3 feature implementation."""

from chemtools.featurizers.calculable import detect_all_features

# Test molecules
test_cases = [
    ("Brc1ccc(Cl)cc1", "4-bromochlorobenzene", {
        "halogen_count": 2,
        "polyhalogenated": True,
    }),
    ("CC(C)(C)c1ccccc1", "tert-butylbenzene", {
        "tert_butyl_present": True,
    }),
    ("CC(C)c1ccccc1", "cumene", {
        "isopropyl_present": True,
    }),
    ("CC(C)(C)OC(=O)Nc1ccc(Br)cc1", "Boc-4-bromoaniline", {
        "boc_present": True,
        "aryl_halide_present": True,
    }),
    ("O=C(NCc1ccccc1)OCc2ccccc2", "Cbz-benzylamine", {
        "cbz_present": True,
    }),
    ("C[Si](C)(C)OCc1ccccc1", "TMS-benzyl ether", {
        "silyl_ether_present": True,
    }),
]

print("="*80)
print("PHASE 3 FEATURE VALIDATION")
print("="*80)

all_passed = True
for smiles, name, expected in test_cases:
    print(f"\n✓ Testing: {name}")
    print(f"  SMILES: {smiles}")
    
    features = detect_all_features(smiles)
    
    for feature, expected_value in expected.items():
        actual_value = features.get(feature)
        status = "✓" if actual_value == expected_value else "✗"
        if actual_value != expected_value:
            all_passed = False
        print(f"  {status} {feature}: {actual_value} (expected {expected_value})")

print("\n" + "="*80)
if all_passed:
    print("✅ ALL VALIDATION TESTS PASSED")
else:
    print("❌ SOME VALIDATION TESTS FAILED")
print("="*80)
