#!/usr/bin/env python3
"""Monkey-patch _knn_impl to add debug output."""

import chemtools.precedent.search as search_module

# Save original function
original_knn_impl = search_module._knn_impl

def debug_knn_impl(family, features, k=50, relax=None):
    """Wrapper with debug output."""
    print(f"\n=== _knn_impl called ===")
    print(f"  family: {repr(family)}")
    print(f"  features: {features}")
    print(f"  k: {k}")
    print(f"  relax keys: {list(relax.keys()) if relax else []}")
    
    result = original_knn_impl(family, features, k, relax)
    
    print(f"  Result: {len(result.get('precedents', []))} precedents")
    if 'error' in result:
        print(f"  Error: {result['error']}")
    print(f"=== end _knn_impl ===\n")
    
    return result

# Monkey-patch
search_module._knn_impl = debug_knn_impl

# Now test
from chemtools import precedent

# Clear caches
from chemtools.precedent.search import _knn_cached
from chemtools.precedent.loader import _load
_knn_cached.cache_clear()
_load.cache_clear()

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

print(f"\nFinal result: {len(pack.get('precedents', []))} precedents")
