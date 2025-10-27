#!/usr/bin/env python3
"""Deep debug: trace through the entire precedent search flow."""

from chemtools.precedent.loader import _load, _load_selective
from chemtools.precedent.search import _candidate_pool

# Clear caches first
_load.cache_clear()

print("=" * 80)
print("STEP 1: Load Suzuki dataset")
print("=" * 80)

rows = _load_selective(['Suzuki'])
print(f"Rows loaded: {len(rows)}")

if rows:
    print(f"\nFirst row:")
    first = rows[0]
    print(f"  reaction_id: {first.get('reaction_id')}")
    print(f"  rxn_type: '{first.get('rxn_type')}'")
    print(f"  reaction_smiles: {first.get('reaction_smiles', 'MISSING')[:60]}...")
    
    # Check rxn_type distribution
    rxn_types = {}
    for row in rows:
        rt = row.get('rxn_type', '(missing)')
        rxn_types[rt] = rxn_types.get(rt, 0) + 1
    
    print(f"\nrxn_type distribution:")
    for rt, count in sorted(rxn_types.items(), key=lambda x: -x[1])[:10]:
        print(f"  '{rt}': {count} reactions")
else:
    print("✗ NO ROWS LOADED!")
    import sys
    sys.exit(1)

print()
print("=" * 80)
print("STEP 2: Build candidate pool")
print("=" * 80)

# Test different family names
for family_test in ['Suzuki', 'suzuki', 'SUZUKI']:
    cands = _candidate_pool(
        rows=rows,
        family_txt=family_test,
        feat={},
        k=20,
        relax={}
    )
    print(f"  Family_txt='{family_test}': {len(cands)} candidates")

print()
print("=" * 80)
print("STEP 3: Check family_text mapping")
print("=" * 80)

from chemtools.precedent.core_utils import _family_text

for name in ['Suzuki', 'suzuki', 'SUZUKI', 'Suzuki-Miyaura', 'suzuki_miyaura']:
    mapped = _family_text(name)
    print(f"  _family_text('{name}') = '{mapped}'")

print()
print("=" * 80)
print("CONCLUSION")
print("=" * 80)

# The most common rxn_type
if rxn_types:
    most_common = max(rxn_types.items(), key=lambda x: x[1])
    print(f"Most common rxn_type: '{most_common[0]}' ({most_common[1]} reactions)")
    
    if most_common[0] != 'Suzuki':
        print(f"\n⚠ FOUND THE BUG!")
        print(f"  Dataset has rxn_type='{most_common[0]}'")
        print(f"  But search is looking for rxn_type='Suzuki'")
        print(f"  → These don't match!")
