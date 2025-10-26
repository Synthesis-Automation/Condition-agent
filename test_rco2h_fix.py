#!/usr/bin/env python
"""Test script to verify rco2h_or_activated_acyl reactant type works correctly."""

from chemtools.analysis import classify_reactant_smiles
from chemtools.analysis.reactants import get_all_reactant_matches

test_cases = [
    ("CC(=O)O", "carboxylic acid"),
    ("CC(=O)Cl", "acyl chloride"),
    ("CC(=O)OC(=O)C", "carboxylic anhydride"),
]

print("Testing rco2h_or_activated_acyl reactant type:")
print("=" * 60)

for smiles, expected_name in test_cases:
    result = classify_reactant_smiles(smiles)
    all_matches = get_all_reactant_matches(smiles)
    
    # Find matches from rco2h_or_activated_acyl
    rco2h_matches = [
        m for m in all_matches 
        if m.member_type in ['RCO2H', 'RCOCl', 'RCO2M', 'Anhydride']
    ]
    
    if rco2h_matches:
        match = rco2h_matches[0]
        print(f"✓ {smiles:20s} -> {match.member_type:10s} ({match.name})")
    else:
        print(f"✗ {smiles:20s} -> No match found!")

print("=" * 60)
print("All tests completed!")
