"""Test Phase 1 critical fixes"""
import sys
sys.path.insert(0, '.')

# Force reload to get latest changes
import importlib
if 'chemtools.featurizers.calculable' in sys.modules:
    importlib.reload(sys.modules['chemtools.featurizers.calculable'])

from chemtools.featurizers.calculable import detect_all_features

print("=" * 70)
print("TESTING PHASE 1 CRITICAL FIXES")
print("=" * 70)

# Test 1: Heteroaryl halides (CRITICAL FIX)
print("\n1. HETEROARYL HALIDE DETECTION (Critical Fix):")
print("-" * 70)
test_heteroaryls = [
    ("4-Bromopyridine", "Brc1ccncc1"),
    ("3-Bromopyridine", "Brc1cccnc1"),
    ("2-Bromopyridine", "Brc1ccccn1"),
    ("2-Bromothiophene", "Brc1cccs1"),
    ("3-Bromofuran", "Brc1ccoc1"),
]

for name, smiles in test_heteroaryls:
    features = detect_all_features(smiles)
    if features:
        het_halide = features.get("heteroaryl_halide_present", False)
        het_present = features.get("heteroaryl_present", False)
        aryl_halide = features.get("aryl_halide_present", False)
        status = "OK  " if het_halide else "FAIL"
        print(f"{status} {name:25s} het_hal: {het_halide}, het: {het_present}, aryl_hal: {aryl_halide}")

# Test 2: Boronic acids vs esters
print("\n2. BORONIC ACID vs ESTER DETECTION:")
print("-" * 70)
test_boronics = [
    ("Phenylboronic acid", "c1ccc(B(O)O)cc1"),
    ("Phenylboronic pinacol ester", "c1ccc(B2OC(C)(C)C(C)(C)O2)cc1"),
]

for name, smiles in test_boronics:
    features = detect_all_features(smiles)
    if features:
        acid = features.get("boronic_acid_present", False)
        ester = features.get("boronic_ester_present", False)
        print(f"{name:30s} acid: {acid}, ester: {ester}")

# Test 3: Primary/Secondary/Tertiary amines
print("\n3. AMINE CLASSIFICATION:")
print("-" * 70)
test_amines = [
    ("Aniline", "c1ccc(N)cc1"),
    ("Benzylamine", "NCc1ccccc1"),
    ("Diethylamine", "CCNCC"),
    ("Triethylamine", "CCN(CC)CC"),
]

for name, smiles in test_amines:
    features = detect_all_features(smiles)
    if features:
        pri = features.get("primary_amine_present", False)
        sec = features.get("secondary_amine_present", False)
        ter = features.get("tertiary_amine_present", False)
        print(f"{name:20s} 1deg: {pri}, 2deg: {sec}, 3deg: {ter}")

# Test 4: Carbonyl features
print("\n4. CARBONYL DETECTION:")
print("-" * 70)
test_carbonyls = [
    ("Acetone", "CC(=O)C"),
    ("Benzaldehyde", "O=Cc1ccccc1"),
    ("Ethyl acetate", "CC(=O)OCC"),
    ("Acetamide", "CC(=O)N"),
]

for name, smiles in test_carbonyls:
    features = detect_all_features(smiles)
    if features:
        ketone = features.get("ketone_present", False)
        aldehyde = features.get("aldehyde_present", False)
        ester = features.get("ester_present", False)
        amide = features.get("amide_present", False)
        print(f"{name:20s} ketone: {ketone}, aldehyde: {aldehyde}, ester: {ester}, amide: {amide}")

# Test 5: Nitrile and Nitro
print("\n5. NITRILE AND NITRO DETECTION:")
print("-" * 70)
test_functional = [
    ("Benzonitrile", "N#Cc1ccccc1"),
    ("Acetonitrile", "CC#N"),
    ("Nitrobenzene", "O=[N+]([O-])c1ccccc1"),
]

for name, smiles in test_functional:
    features = detect_all_features(smiles)
    if features:
        nitrile = features.get("nitrile_present", False)
        aryl_nitrile = features.get("aryl_nitrile_present", False)
        nitro = features.get("nitro_present", False)
        aryl_nitro = features.get("aryl_nitro_present", False)
        print(f"{name:20s} nitrile: {nitrile}, aryl_N: {aryl_nitrile}, nitro: {nitro}, aryl_NO2: {aryl_nitro}")

print("\n" + "=" * 70)
print("PHASE 1 TESTING COMPLETE")
print("=" * 70)
