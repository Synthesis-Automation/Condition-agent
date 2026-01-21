#!/usr/bin/env python3
"""
Remove redundant 'name' field from organic_compounds.v1.3.json

Since we simplified the system so that name == id, the 'name' field is redundant.
This script removes it to improve maintainability.
"""

import json
from pathlib import Path

# Get the absolute path
taxonomy_dir = Path(__file__).parent / "data"
compounds_file = taxonomy_dir / "organic_compounds.v1.3.json"

# Load
print(f"Loading {compounds_file.name}...")
with open(compounds_file, 'r', encoding='utf-8') as f:
    data = json.load(f)

compounds = data.get('compounds', [])
print(f"Found {len(compounds)} compounds")

# Remove name field and verify it matches id
removed = 0
mismatches = []

for compound in compounds:
    compound_id = compound.get('id')
    compound_name = compound.get('name')
    
    if 'name' in compound:
        # Verify name matches id before removing
        if compound_name != compound_id:
            mismatches.append({
                'id': compound_id,
                'name': compound_name
            })
        
        del compound['name']
        removed += 1

# Report
print(f"\n{'='*70}")
print(f"Removed 'name' field from {removed} compounds")

if mismatches:
    print(f"\n⚠ WARNING: Found {len(mismatches)} compounds where name != id:")
    for m in mismatches:
        print(f"  - {m['id']}: name was '{m['name']}'")
    print("\nThese should have been fixed during simplification!")
else:
    print("✓ All names matched their IDs (as expected)")

# Save
with open(compounds_file, 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=2, ensure_ascii=False)

print(f"\n✓ Updated {compounds_file.name}")
print(f"\nBenefits:")
print("  - Reduced file size")
print("  - Single source of truth (ID only)")
print("  - No chance of name/ID mismatch")
print("  - Easier to maintain")
print("  - Cleaner JSON structure")
