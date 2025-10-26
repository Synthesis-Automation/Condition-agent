#!/usr/bin/env python
"""Test rxn_insight detection on Phase 2 reaction types."""

from chemtools.reaction_type_detector import detect_reaction_type
from chemtools.router import detect_family_from_reaction

test_reactions = [
    # Esterification
    ("CC(=O)O.CCO>>CC(=O)OCC", "Esterification (Fischer)"),
    ("c1ccc(C(=O)O)cc1.OCC>>c1ccc(C(=O)OCC)cc1", "Esterification (benzoic acid)"),
    
    # Grignard addition
    ("c1ccccc1C(=O)C.C[Mg]Br>>c1ccccc1C(C)(C)O", "Grignard addition to ketone"),
    
    # Hydroboration
    ("C=Cc1ccccc1.B2H6>>CCc1ccccc1", "Hydroboration of styrene"),
    
    # Nitrile formation (SN2 with CN-)
    ("BrCCCC.[CN-]>>N#CCCCC", "Nitrile formation (SN2)"),
    
    # Finkelstein (halide exchange)
    ("BrCCCC.[I-]>>ICCCC", "Finkelstein (Br→I)"),
    
    # Williamson ether
    ("BrCCCC.CCO>>CCOCCCC", "Williamson ether synthesis"),
    
    # Claisen condensation
    ("CC(=O)OCC.CC(=O)OCC>>CC(=O)CC(=O)OCC", "Claisen condensation"),
    
    # Michael addition
    ("C=CC(=O)C.CC(=O)CC(=O)OCC>>CCC(C(=O)C)C(C(=O)OCC)C(=O)C", "Michael addition"),
]

print("=" * 80)
print("PHASE 2 REACTION TYPE DETECTION TEST")
print("=" * 80)

for rxn_smiles, description in test_reactions:
    print(f"\n{description}")
    print(f"SMILES: {rxn_smiles}")
    
    # Test rule-based
    rule_result = detect_family_from_reaction(rxn_smiles, use_rxn_insight=False)
    print(f"  Rule-based: {rule_result.get('family')} (conf: {rule_result.get('confidence'):.2f})")
    
    # Test with rxn_insight
    full_result = detect_family_from_reaction(rxn_smiles, use_rxn_insight=True)
    auto = full_result.get('auto')
    if auto and auto.get('success'):
        print(f"  rxn_insight: {auto.get('rxn_class')} → {auto.get('mapped_family')} (conf: {auto.get('confidence') or 'N/A'})")
    else:
        print(f"  rxn_insight: Not detected")
    
    print(f"  FINAL: {full_result.get('family')}")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("Current detection method: Rule-based (7 patterns) + rxn_insight (ML, optional)")
print("Rule-based covers: C-N, C-O, C-S, Suzuki, Sonogashira, Amide coupling")
print("rxn_insight: Installed but may not accurately detect all Phase 2 types")
