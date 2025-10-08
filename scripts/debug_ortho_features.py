#!/usr/bin/env python3
"""Debug script to check ortho feature extraction."""

import sys
sys.path.insert(0, ".")

from chemtools.scdb_matcher import matcher, loader

# Test reaction with 2-bromotoluene (1 ortho substituent)
rxn_smiles = "Brc1ccccc1C.c1ccc(B(O)O)cc1>>Cc1ccccc1-c1ccccc1"

# Load database
db = loader.load_db("data/conditionDB/suzuki_db.json")

# Parse reactants and get features
reactants = matcher._parse_reactant_smiles(rxn_smiles)
print(f"Reactants: {reactants}")

normalization = matcher.essential_core_normalize(reactants)
print(f"\nNormalized reactants: {normalization.kept_reactants}")
print(f"Masked mols: {len(normalization.masked_mols)} molecules")

set_features, numeric_features = matcher._compile_feature_sets(normalization)

print("\n" + "="*70)
print("EXTRACTED FEATURES")
print("="*70)

print("\nSet Features:")
for key, values in sorted(set_features.items()):
    print(f"  {key}: {values}")

print("\nNumeric Features:")
for key, value in sorted(numeric_features.items()):
    print(f"  {key}: {value}")

print("\n" + "="*70)
print("CHECKING XPHOS ENTRY REQUIREMENTS")
print("="*70)

# Find XPhos entry
xphos_entry = None
for entry in db.entries:
    if entry.id == "SCDB-SUZ-ARBR-ORTHO-XPhos":
        xphos_entry = entry
        break

if xphos_entry:
    print(f"\nEntry ID: {xphos_entry.id}")
    print(f"Feature Requirements: {xphos_entry.feature_requirements}")
    
    print("\nChecking requirements:")
    if xphos_entry.feature_requirements:
        for key, expected in xphos_entry.feature_requirements.items():
            if key in numeric_features:
                actual = numeric_features[key]
                print(f"  {key}: actual={actual}, expected={expected}")
                result = matcher._check_numeric_requirement(actual, expected)
                print(f"    ‚Ü?Satisfied: {result}")
            elif key in set_features:
                actual = set_features[key]
                print(f"  {key}: actual={actual}, expected={expected}")
                result = matcher._satisfies_set_requirement(actual, expected)
                print(f"    ‚Ü?Satisfied: {result}")
            else:
                print(f"  {key}: NOT FOUND in features")
                print(f"    ‚Ü?Satisfied: False")
    
    applicable = matcher._entry_applicable(xphos_entry, set_features, numeric_features)
    print(f"\nEntry is applicable: {applicable}")
else:
    print("\n‚ù?SCDB-SUZ-ARBR-ORTHO-XPhos not found in database!")

print("\n" + "="*70)
