#!/usr/bin/env python3
"""Debug why precedent search returns 0 results."""

import sys
import json

# Check what's in the Suzuki dataset
print("=" * 80)
print("INSPECTING SUZUKI DATASET")
print("=" * 80)

# Read first few lines of Suzuki.jsonl to see the structure
with open('data/reaction_dataset/Suzuki.jsonl', 'r', encoding='utf-8') as f:
    lines = f.readlines()[:5]  # First 5 reactions
    
print(f"Total reactions in file: {len(open('data/reaction_dataset/Suzuki.jsonl', 'r', encoding='utf-8').readlines())}")
print(f"\nFirst 5 reactions:")
print()

for i, line in enumerate(lines, 1):
    rxn = json.loads(line)
    print(f"Reaction {i}:")
    print(f"  rxn_type: '{rxn.get('rxn_type', 'MISSING')}'")
    print(f"  reaction_smiles: {rxn.get('reaction_smiles', 'MISSING')[:60]}...")
    if 'features' in rxn:
        print(f"  features.bin: {rxn['features'].get('bin', 'MISSING')}")
    print()

# Now check what precedent._load() returns
print("=" * 80)
print("CHECKING precedent._load()")
print("=" * 80)

from chemtools.precedent.search import _load, _load_selective

# Load all datasets
print("Loading Suzuki data via _load_selective(['Suzuki'])...")
rows = _load_selective(['Suzuki'])
print(f"Rows loaded: {len(rows)}")

if rows:
    print(f"\nFirst row structure:")
    first_row = rows[0]
    print(f"  rxn_type: '{first_row.get('rxn_type', 'MISSING')}'")
    print(f"  reaction_smiles: {first_row.get('reaction_smiles', 'MISSING')[:60]}...")
    if 'features' in first_row:
        print(f"  features.bin: {first_row['features'].get('bin', 'MISSING')}")
    
    # Check rxn_type distribution
    rxn_types = {}
    for row in rows:
        rt = row.get('rxn_type', '(missing)')
        rxn_types[rt] = rxn_types.get(rt, 0) + 1
    
    print(f"\nrxn_type distribution:")
    for rt, count in sorted(rxn_types.items(), key=lambda x: -x[1])[:10]:
        print(f"  '{rt}': {count} reactions")
else:
    print("NO ROWS LOADED!")

print()
print("=" * 80)
print("TESTING CANDIDATE POOL")
print("=" * 80)

from chemtools.precedent.search import _candidate_pool

# Test what _candidate_pool returns
cands = _candidate_pool(
    rows=rows,
    family='Suzuki',
    features={},
    k=20,
    relax={}
)

print(f"Candidates returned: {len(cands)}")

if not cands:
    print("\nDEBUG: Why no candidates?")
    print("  Checking family filter...")
    
    # Try different family names
    for family_test in ['Suzuki', 'suzuki', 'SUZUKI', 'Suzuki-Miyaura', 'suzuki_miyaura']:
        filtered = [r for r in rows if r.get('rxn_type') == family_test]
        if filtered:
            print(f"  '{family_test}': {len(filtered)} matches ✓")
        else:
            print(f"  '{family_test}': 0 matches")
    
    # Check actual rxn_type values
    print(f"\n  Actual rxn_type values in dataset:")
    unique_types = set(r.get('rxn_type', '(missing)') for r in rows[:100])
    for rt in sorted(unique_types):
        print(f"    '{rt}'")

print()
print("=" * 80)
print("CONCLUSION")
print("=" * 80)
if cands:
    print("✓ Candidate pool building is working")
else:
    print("✗ Candidate pool is EMPTY - this is the bug!")
    print("  → Check rxn_type field in Suzuki.jsonl")
    print("  → Family name mismatch?")
