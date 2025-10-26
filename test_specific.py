#!/usr/bin/env python
"""Test specific reactions to see what's being detected."""

import sys
sys.path.insert(0, '.')

from chemtools.router import detect_family_from_reaction

test_cases = [
    # Heck coupling
    ("Brc1ccccc1.C=Cc1ccccc1>>c1ccccc1C=Cc2ccccc2", "Heck"),
    # Kumada coupling  
    ("Brc1ccccc1.C[Mg]Br>>Cc1ccccc1", "Kumada"),
    # Negishi coupling
    ("Brc1ccccc1.C[Zn]Cl>>Cc1ccccc1", "Negishi"),
    # Hydrogenation (should only override if no better match)
    ("C=Cc1ccccc1.[H][H]>>CCc1ccccc1", "Hydrogenation"),
    # Williamson ether (alkyl halide + alcohol)
    ("BrCCCC.CCO>>CCOCCCC", "Williamson"),
]

for rxn_smiles, expected in test_cases:
    result = detect_family_from_reaction(rxn_smiles, use_rxn_insight=False)
    family = result.get('family', 'UNKNOWN')
    conf = result.get('confidence', 0.0)
    print(f"{expected:20s} → {family:25s} (conf: {conf:.2f})")
