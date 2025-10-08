#!/usr/bin/env python3
"""
Debug SMARTS Matching
=====================
Test if SMARTS patterns actually match the test molecules.
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from rdkit import Chem

# Test case: SCDB-SUZ-ARBRI-GENERAL-PPh3
rxn_smiles = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
smarts_pattern = "[c:1]-[I,Br:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]"  # FIXED!

print("Testing SMARTS pattern matching")
print("=" * 70)
print(f"Reaction SMILES: {rxn_smiles}")
print(f"SMARTS pattern:  {smarts_pattern}")
print()

# Split reaction
reactants, products = rxn_smiles.split(">>")
reactant_smiles = [s.strip() for s in reactants.split(".")]

print(f"Reactants: {reactant_smiles}")
print()

# Split SMARTS
smarts_reactants = smarts_pattern.split(">>")[0]
smarts_patterns = [s.strip() for s in smarts_reactants.split(".")]

print(f"SMARTS patterns: {smarts_patterns}")
print()

# Test each reactant against each pattern
for i, reactant in enumerate(reactant_smiles):
    mol = Chem.MolFromSmiles(reactant)
    if mol is None:
        print(f"Reactant {i}: INVALID SMILES '{reactant}'")
        continue
    
    print(f"Reactant {i}: {reactant}")
    for j, pattern_str in enumerate(smarts_patterns):
        pattern = Chem.MolFromSmarts(pattern_str)
        if pattern is None:
            print(f"  Pattern {j}: INVALID SMARTS '{pattern_str}'")
            continue
        
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            print(f"  Pattern {j}: MATCH ({len(matches)} matches) - {pattern_str}")
        else:
            print(f"  Pattern {j}: NO MATCH - {pattern_str}")
    print()

# Try the boronic acid pattern specifically
print("Testing boronic acid pattern variations:")
print("-" * 70)
boronic_acid = "OB(O)c1ccccc1"
mol = Chem.MolFromSmiles(boronic_acid)

test_patterns = [
    "[B:3]([O])[O]-[c:4]",  # Current pattern
    "[B:3](-[O])(-[O])-[c:4]",  # Explicit
    "c-[B]([O])[O]",  # No atom mapping
    "c-B(O)O",  # Simple
]

for pattern_str in test_patterns:
    pattern = Chem.MolFromSmarts(pattern_str)
    if pattern is None:
        print(f"  INVALID: {pattern_str}")
        continue
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        print(f"  MATCH: {pattern_str} ({len(matches)} matches)")
    else:
        print(f"  NO MATCH: {pattern_str}")
