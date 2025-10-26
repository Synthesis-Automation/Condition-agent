#!/usr/bin/env python
"""Debug why reactions are still UNKNOWN after Phase 2."""

import sys
sys.path.insert(0, '.')

from chemtools.analysis import analyze_reaction
from chemtools.router import detect_family_from_reaction

# Test some specific UNKNOWN reactions
test_cases = [
    ("CC(=O)O.CCO>>CC(=O)OCC", "Esterification"),
    ("c1ccccc1C=C.B2H6>>CCc1ccccc1", "Hydroboration"),
    ("BrCCCC.[C-]#N>>N#CCCCC", "Nitrile formation"),
    ("BrCCCC.[I-]>>ICCCC", "Finkelstein"),
    ("c1ccccc1Br.C[Mg]Br>>c1ccccc1C", "Grignard (should be detected)"),
    ("C=Cc1ccccc1.[H][H]>>CCc1ccccc1", "Hydrogenation - H2"),
    ("c1ccccc1C(=O)C.B2H6>>c1ccccc1C(C)O", "Carbonyl reduction - BH3"),
    ("CCO.CrO3>>CC=O", "Alcohol oxidation"),
]

print("=" * 100)
print("DEBUGGING PHASE 2 DETECTION")
print("=" * 100)

for rxn_smiles, expected in test_cases:
    print(f"\n{expected}")
    print(f"  SMILES: {rxn_smiles}")
    
    result = analyze_reaction(rxn_smiles, use_rxn_insight=False)
    detection = result.get('family', {}).get('detected', {})
    
    family = detection.get('family', 'UNKNOWN')
    hits = detection.get('hits', {})
    
    print(f"  Detected: {family}")
    print(f"  Hits: {', '.join([k for k,v in hits.items() if v])}")
    
    # Check agents
    agents = result.get('agents', [])
    if agents:
        print(f"  Agents: {agents}")
    else:
        print(f"  Agents: (none detected)")
    
    # Check reactants
    reactants = result.get('reactants', [])
    print(f"  Reactants ({len(reactants)}):") 
    for i, r in enumerate(reactants):
        norm = r.get('normalized', {})
        best = r.get('taxonomy', {}).get('best_match', {})
        print(f"    [{i+1}] {norm.get('smiles_norm', '?')} → {best.get('name', '?')}")

print("\n" + "=" * 100)
