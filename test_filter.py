#!/usr/bin/env python3
"""Test if filter_by_reagent_database is the problem."""

from chemtools import precedent
from chemtools.precedent.search import _knn_cached
from chemtools.precedent.loader import _load

# Clear caches
_knn_cached.cache_clear()
_load.cache_clear()

reaction_smiles = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"

print("=" * 80)
print("TEST 1: WITH filter_by_reagent_database=True (default)")
print("=" * 80)
pack1 = precedent.knn(
    family='Suzuki',
    features={},
    k=20,
    relax={
        'use_drfp': True,
        'reaction_smiles': reaction_smiles,
        'filter_by_reagent_database': True
    }
)
print(f"Precedents: {len(pack1.get('precedents', []))}")
print()

print("=" * 80)
print("TEST 2: WITHOUT filter_by_reagent_database (set to False)")
print("=" * 80)
_knn_cached.cache_clear()  # Clear cache to force recomputation
pack2 = precedent.knn(
    family='Suzuki',
    features={},
    k=20,
    relax={
        'use_drfp': True,
        'reaction_smiles': reaction_smiles,
        'filter_by_reagent_database': False
    }
)
print(f"Precedents: {len(pack2.get('precedents', []))}")

if pack2.get('precedents'):
    print(f"\n✓ SUCCESS! The reagent database filter was removing all precedents!")
    print(f"\nFirst 3 precedents:")
    for i, prec in enumerate(pack2['precedents'][:3], 1):
        print(f"\n{i}. {prec.get('reaction_id')}")
        print(f"   Catalyst: {prec.get('catalyst_name', 'N/A')}")
        print(f"   Base: {prec.get('base_name', 'N/A')}")
        print(f"   Solvent: {prec.get('solvent_name', 'N/A')}")
else:
    print(f"\n✗ Still 0 precedents - the filter is NOT the problem")
