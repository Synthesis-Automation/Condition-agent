#!/usr/bin/env python
"""Debug how agents are parsed in reactions."""

import sys
sys.path.insert(0, '.')

from chemtools.smiles import normalize_reaction

# Test reactions with known reagents
test_cases = [
    ("C=Cc1ccccc1.[H][H]>>CCc1ccccc1", "Hydrogenation with H2"),
    ("c1ccccc1C(=O)C>>c1ccccc1C(C)O", "Reduction (should have NaBH4)"),
    ("CCO>>CC=O", "Oxidation (should have oxidant)"),
    ("C=Cc1ccccc1.B2H6>>CCc1ccccc1", "Hydroboration with B2H6"),
]

print("=" * 100)
print("AGENT PARSING DEBUG")
print("=" * 100)

for rxn_smiles, desc in test_cases:
    print(f"\n{desc}")
    print(f"  SMILES: {rxn_smiles}")
    
    norm = normalize_reaction(rxn_smiles)
    
    print(f"  Normalized: {norm.get('normalized', '?')}")
    
    reactants = norm.get('reactants', [])
    print(f"  Reactants ({len(reactants)}): {[r.get('smiles_norm', '?') for r in reactants]}")
    
    agents = norm.get('agents', [])
    if agents:
        print(f"  Agents ({len(agents)}): {[a.get('smiles_norm', '?') for a in agents]}")
    else:
        print(f"  Agents: (empty list)")
    
    products = norm.get('products', [])
    print(f"  Products ({len(products)}): {[p.get('smiles_norm', '?') for p in products]}")

print("\n" + "=" * 100)
print("ISSUE: Reagents like H2, NaBH4 appear in reactants, not agents!")
print("Solution: Check all reactants for common reagents and flag them")
print("=" * 100)
