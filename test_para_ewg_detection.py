"""
Verify that para_EWG detection works correctly for aryl halides.
"""

from chemtools.featurizers.molecular import featurize

test_cases = [
    {
        "name": "p-Anisidine (EDG -OMe)",
        "electrophile": "O=Cc1ccccc1",
        "nucleophile": "Nc1ccc(OC)cc1",
        "expected_para_EWG": False,
        "reason": "Benzaldehyde has no leaving group"
    },
    {
        "name": "p-Nitro-chlorobenzene (EWG -NO2)",
        "electrophile": "Clc1ccc([N+](=O)[O-])cc1",
        "nucleophile": "Nc1ccccc1",
        "expected_para_EWG": True,
        "reason": "Aryl chloride with para-nitro group"
    },
    {
        "name": "p-Chlorobenzonitrile (EWG -CN)",
        "electrophile": "Clc1ccc(C#N)cc1",
        "nucleophile": "Nc1ccccc1",
        "expected_para_EWG": True,
        "reason": "Aryl chloride with para-cyano group"
    },
    {
        "name": "Simple chlorobenzene (no EWG)",
        "electrophile": "Clc1ccccc1",
        "nucleophile": "Nc1ccccc1",
        "expected_para_EWG": False,
        "reason": "Aryl chloride without EWG"
    },
    {
        "name": "p-Trifluoromethyl-bromobenzene (EWG -CF3)",
        "electrophile": "Brc1ccc(C(F)(F)F)cc1",
        "nucleophile": "Nc1ccccc1",
        "expected_para_EWG": True,
        "reason": "Aryl bromide with para-CF3 group"
    }
]

print("=" * 80)
print("Testing para_EWG Detection for Aryl Halides")
print("=" * 80)
print()

all_passed = True

for test in test_cases:
    print(f"Test: {test['name']}")
    print(f"  Electrophile: {test['electrophile']}")
    print(f"  Nucleophile: {test['nucleophile']}")
    
    features = featurize(test['electrophile'], test['nucleophile'])
    detected_para_EWG = features.get('para_EWG', False)
    detected_LG = features.get('LG', 'UNK')
    detected_class = features.get('elec_class', 'UNK')
    
    print(f"  LG: {detected_LG}, elec_class: {detected_class}")
    print(f"  para_EWG: {detected_para_EWG} (expected: {test['expected_para_EWG']})")
    
    if detected_para_EWG == test['expected_para_EWG']:
        print(f"  ✅ PASS - {test['reason']}")
    else:
        print(f"  ❌ FAIL - {test['reason']}")
        all_passed = False
    print()

print("=" * 80)
if all_passed:
    print("✅ ALL TESTS PASSED")
else:
    print("❌ SOME TESTS FAILED")
print("=" * 80)
