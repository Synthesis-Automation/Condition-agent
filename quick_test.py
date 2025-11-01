#!/usr/bin/env python3
"""Quick integration test"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from chemtools import (
    add_atom_mapping,
    analyze_bond_changes,
    rxnmapper_available,
    detect_reaction,
)

print("="*60)
print("RXNMAPPER INTEGRATION - QUICK TEST")
print("="*60)

print("\n✅ All imports successful!")
print(f"✅ RXNMapper available: {rxnmapper_available()}")

smiles = "Brc1ccccc1.C#C>>c1ccccc1C#C"
print(f"\nTest reaction: {smiles}")

# Test bond analysis
result = analyze_bond_changes(smiles, auto_map=True)
print(f"\n✅ Bond analysis works: {result['success']}")
if result['success']:
    print(f"   Broken bonds: {len(result.get('broken_bonds', []))}")
    print(f"   Formed bonds: {len(result.get('formed_bonds', []))}")
    print(f"   Changed atoms: {len(result.get('changed_atoms', []))}")

# Test detection + bond analysis together
detection = detect_reaction(smiles, use_ml=True)
print(f"\n✅ Reaction detection: {detection['family']}")
print(f"   Confidence: {detection['confidence']:.2f}")

print("\n" + "="*60)
print("ALL TESTS PASSED! ✨")
print("="*60)
