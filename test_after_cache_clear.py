#!/usr/bin/env python3
"""Clear the precedent search cache and retest."""

from chemtools.precedent.search import _knn_cached

# Clear the LRU cache
print("Clearing precedent search cache...")
_knn_cached.cache_clear()
print("✓ Cache cleared")
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
    print(f"\nFirst precedent:")
    first = pack['precedents'][0]
    print(f"  reaction_id: {first.get('reaction_id')}")
    print(f"  catalyst: {first.get('catalyst_name')}")
    print(f"  base: {first.get('base_name')}")
    print(f"  solvent: {first.get('solvent_name')}")
else:
    print(f"\n✗ Still getting 0 precedents")
    if 'error' in pack:
        print(f"  Error: {pack['error']}")
