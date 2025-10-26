#!/usr/bin/env python
"""Debug why Heck/Kumada/Negishi aren't detected."""

import sys
sys.path.insert(0, '.')

from chemtools.analysis import analyze_reaction

test_cases = [
    ("Brc1ccccc1.C=C>>C=Cc1ccccc1", "Heck - Simple styrene"),
    ("Brc1ccccc1.C[Mg]Br>>Cc1ccccc1", "Kumada - Ph-Toluene"),  
    ("Brc1ccccc1.C[Zn]Br>>Cc1ccccc1", "Negishi - Ph-Toluene"),
]

for rxn_smiles, expected in test_cases:
    print(f"\n{expected}")
    print(f"  SMILES: {rxn_smiles}")
    
    result = analyze_reaction(rxn_smiles, use_rxn_insight=False)
    detection = result.get('family', {}).get('detected', {})
    
    family = detection.get('family', 'UNKNOWN')
    hits = detection.get('hits', {})
    
    print(f"  Detected: {family}")
    print(f"  Relevant hits:")
    for k in ['aryl_halide', 'alkene', 'grignard', 'organozinc']:
        if k in hits:
            print(f"    {k}: {hits[k]}")
