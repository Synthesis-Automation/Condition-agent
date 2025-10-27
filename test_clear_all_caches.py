#!/usr/bin/env python3
"""Clear ALL caches and retest."""

from chemtools.precedent.search import _knn_cached
from chemtools.precedent.loader import _load

# Clear both caches
print("Clearing all precedent caches...")
_knn_cached.cache_clear()
_load.cache_clear()
print("✓ All caches cleared")
print()

# Now re-run the test
from chemtools import precedent

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

if pack.get('precedents'):
    print(f"\n✓ SUCCESS! Precedent search is working now!")
    print(f"\nFirst 3 precedents:")
    for i, prec in enumerate(pack['precedents'][:3], 1):
        print(f"\n{i}. {prec.get('reaction_id')}")
        print(f"   Catalyst: {prec.get('catalyst_name', 'N/A')}")
        print(f"   Base: {prec.get('base_name', 'N/A')}")
        print(f"   Solvent: {prec.get('solvent_name', 'N/A')}")
        print(f"   Yield: {prec.get('yield_value', 'N/A')}%")
else:
    print(f"\n✗ Still getting 0 precedents")
    if 'error' in pack:
        print(f"  Error: {pack['error']}")
