#!/usr/bin/env python3
"""Test precedent.knn() directly to isolate the issue."""

from chemtools import precedent

# Test 1: Suzuki with explicit DRFP
print("=" * 80)
print("TEST 1: Suzuki with DRFP")
print("=" * 80)

reaction_smiles = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
pack = precedent.knn(
    family='Suzuki',
    features={},
    k=20,
    relax={
        'use_drfp': True,
        'reaction_smiles': reaction_smiles
    }
)

print(f"Reaction: {reaction_smiles}")
print(f"Precedents returned: {len(pack.get('precedents', []))}")
print(f"Support: {pack.get('support', 0)}")
if 'error' in pack:
    print(f"Error: {pack['error']}")
print()

# Test 2: Suzuki WITHOUT DRFP (feature-based only)
print("=" * 80)
print("TEST 2: Suzuki WITHOUT DRFP (features only)")
print("=" * 80)

pack2 = precedent.knn(
    family='Suzuki',
    features={},
    k=20,
    relax={
        'use_drfp': False,
        'reaction_smiles': reaction_smiles
    }
)

print(f"Reaction: {reaction_smiles}")
print(f"Precedents returned: {len(pack2.get('precedents', []))}")
print(f"Support: {pack2.get('support', 0)}")
if 'error' in pack2:
    print(f"Error: {pack2['error']}")
print()

# Test 3: Suzuki with DRFP but NO reaction_smiles in relax
print("=" * 80)
print("TEST 3: Suzuki with DRFP but NO reaction_smiles")
print("=" * 80)

pack3 = precedent.knn(
    family='Suzuki',
    features={},
    k=20,
    relax={
        'use_drfp': True
        # NO reaction_smiles!
    }
)

print(f"Precedents returned: {len(pack3.get('precedents', []))}")
print(f"Support: {pack3.get('support', 0)}")
if 'error' in pack3:
    print(f"Error: {pack3['error']}")
print()

# Test 4: Compare with how recommender calls it
print("=" * 80)
print("TEST 4: Call pattern from recommender.py")
print("=" * 80)

relax_from_recommender = {
    'reaction_smiles': reaction_smiles,
    'use_drfp': True,
    'precompute_drfp': False,
    'selective_loading': True
}

pack4 = precedent.knn(
    family='Suzuki',
    features={},
    k=20,
    relax=relax_from_recommender
)

print(f"Reaction: {reaction_smiles}")
print(f"Precedents returned: {len(pack4.get('precedents', []))}")
print(f"Support: {pack4.get('support', 0)}")
if 'error' in pack4:
    print(f"Error: {pack4['error']}")
print()

print("=" * 80)
print("CONCLUSION")
print("=" * 80)
if len(pack.get('precedents', [])) > 0:
    print("✓ Precedent search IS working when called directly")
    print("  → The bug must be in how recommender.py calls it")
else:
    print("✗ Precedent search is broken even when called directly")
    print("  → The bug is in precedent.knn() itself")
