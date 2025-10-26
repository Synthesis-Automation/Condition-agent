#!/usr/bin/env python
"""Debug specific reactions to see why they're UNKNOWN."""

import sys
sys.path.insert(0, '.')

from chemtools.router import detect_family_from_reaction
from chemtools.analysis import normalize_reaction

# Test reactions from sample_reactions.py
test_cases = [
    # Heck (should be detected)
    ("C=Cc1ccccc1.Ic1ccccc1>>C=Cc1ccc(cc1)c1ccccc1", "Heck - Simple styrene"),
    
    # Kumada (should be detected)
    ("Brc1ccccc1.[Mg]Br.c1ccccc1>>c1ccc(cc1)c1ccccc1", "Kumada - Ph-Ph"),
    
    # Negishi (should be detected)
    ("Brc1ccccc1.[Zn].c1ccccc1>>c1ccc(cc1)c1ccccc1", "Negishi - Ph-Ph"),
    
    # Hydrogenation (should be detected, was failing before fix)
    ("c1ccc(C=O)cc1.[H][H]>>c1ccc(CO)cc1", "Hydrogenation - Benzyl alcohol"),
    
    # Williamson ether (should be detected)
    ("Brc1ccccc1.c1ccc(O)cc1>>c1ccc(Oc2ccccc2)cc1", "Williamson ether"),
]

print("Testing specific reactions:\n")
print("=" * 100)

for smiles, desc in test_cases:
    print(f"\n{desc}")
    print(f"SMILES: {smiles}")
    
    # First check normalization
    norm = normalize_reaction(smiles)
    print(f"Normalized: {norm.get('normalized', 'ERROR')}")
    print(f"Reactants: {len(norm.get('reactants', []))} molecules")
    
    # Get reactant SMILES
    reactant_smiles = []
    for r in norm.get('reactants', []):
        s = r.get('smiles_norm') or r.get('largest_smiles') or r.get('input', '')
        if s:
            reactant_smiles.append(s)
    
    print(f"Reactant SMILES: {'.'.join(reactant_smiles)}")
    
    # Detect family - direct call to detect_family()
    from chemtools.router import detect_family
    base_result = detect_family(reactant_smiles)
    print(f"Base detect_family(): {base_result}")
    
    # Detect family - call through detect_family_from_reaction()
    result = detect_family_from_reaction(smiles, use_rxn_insight=False)
    family = result.get('canonical_id', 'UNKNOWN')
    family_direct = result.get('family', 'UNKNOWN')
    confidence = result.get('confidence', 0.0)
    
    print(f"detect_family_from_reaction() result: {result}")
    print(f"FINAL DETECTED: family={family}, family_direct={family_direct}, confidence={confidence:.2f}")
    print("=" * 100)
