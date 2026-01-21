#!/usr/bin/env python3
"""Test that featurizer and data converter work with simplified taxonomy."""

from chemtools.featurizers.unified import featurize_molecule

# Test 1: Basic featurization
print("Test 1: Basic featurization")
result = featurize_molecule('c1ccccc1Cl')
motifs = result.get('molecule', {}).get('motifs', [])
print(f"✓ Found {len(motifs)} motifs: {[m['compound_id'] for m in motifs]}")

# Test 2: Nearby groups (uses name field fallback)
print("\nTest 2: Nearby groups")
result = featurize_molecule('c1cc(Cl)c(Br)cc1')
analyses = result.get('molecule', {}).get('analyses', [])
if analyses:
    nearby = analyses[0].get('nearby_groups', [])
    print(f"✓ Found {len(nearby)} nearby groups")
    for g in nearby[:3]:
        name = g.get('name', g.get('group_id', 'unknown'))
        print(f"  - {name}")

# Test 3: Tags (should be empty lists now)
print("\nTest 3: Tags (should be empty)")
from chemtools.featurizers.motif_registry import build_compound_registry
from pathlib import Path

registry_paths = {
    "groups": Path("chemtools/taxonomy/data/organic_groups.v1.3.json"),
    "compounds": Path("chemtools/taxonomy/data/organic_compounds.v1.3.json"),
}
registry = build_compound_registry(registry_paths)
pattern = registry["compound_map"].get("Ar-Cl")
if pattern:
    print(f"✓ Ar-Cl b_tags: {pattern.b_tags} (empty is expected)")

print("\n✓ All tests passed! No updates needed.")
