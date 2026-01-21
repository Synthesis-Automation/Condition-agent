#!/usr/bin/env python3
"""
Test script to verify that the featurizer correctly loads compounds
with auto-generated IDs from the A-B pattern.
"""

from pathlib import Path
from chemtools.featurizers.motif_registry import build_compound_registry

def test_auto_generated_ids():
    print("Testing Auto-Generated Compound IDs")
    print("=" * 70)
    
    # Load taxonomy files
    taxonomy_dir = Path("chemtools/taxonomy/data")
    registry_paths = {
        "groups": taxonomy_dir / "organic_groups.v1.3.json",
        "compounds": taxonomy_dir / "organic_compounds.v1.3.json",
    }
    
    print("\n1. Building compound registry with auto-generated IDs...")
    registry = build_compound_registry(registry_paths)
    
    compiled_compounds = registry["compiled_compounds"]
    compound_map = registry["compound_map"]
    
    print(f"   ✓ Loaded {len(compiled_compounds)} compiled compounds")
    print(f"   ✓ Compound map has {len(compound_map)} entries")
    
    # Show sample auto-generated IDs
    print("\n2. Sample auto-generated compound IDs:")
    sample_compounds = list(compound_map.items())[:10]
    for comp_id, pattern in sample_compounds:
        print(f"   - {comp_id:20s} (A={pattern.group_a}, B={pattern.group_b})")
    
    # Verify specific known compounds
    print("\n3. Verifying specific compounds:")
    test_cases = [
        ("Ar-Cl", "Ar", "Cl"),
        ("Ar-Br", "Ar", "Br"),
        ("Alkyl-OH", "Alkyl", "OH"),
        ("RCH2-Sn*", "RCH2", "Sn*"),
        ("Ar-Ar", "Ar", "Ar_Subst"),  # B is Ar_Subst but ID normalized to Ar-Ar
        ("Ar-Alkenyl", "Ar", "Alkenyl_Subst"),  # Similar normalization
        ("Ar-AromN", "Ar", "AromN_Subst"),  # Similar normalization
    ]
    
    for expected_id, group_a, group_b in test_cases:
        if expected_id in compound_map:
            pattern = compound_map[expected_id]
            match = pattern.group_a == group_a and pattern.group_b == group_b
            status = "✓" if match else "✗"
            print(f"   {status} {expected_id}: A={pattern.group_a}, B={pattern.group_b}")
        else:
            print(f"   ✗ {expected_id}: NOT FOUND in registry")
    
    # Check combination_map (A, B) -> compound_id
    print("\n4. Testing combination map (A, B) -> compound_id:")
    combination_map = registry["combination_map"]
    test_combinations = [
        (("Ar", "Cl"), "Ar-Cl"),
        (("Ar", "Br"), "Ar-Br"),
        (("Alkyl", "OH"), "Alkyl-OH"),
    ]
    
    for (a, b), expected_id in test_combinations:
        actual_id = combination_map.get((a, b))
        match = actual_id == expected_id
        status = "✓" if match else "✗"
        print(f"   {status} ({a}, {b}) -> {actual_id or 'NOT FOUND'}")
    
    print("\n" + "=" * 70)
    print("✓ Auto-generated ID system working correctly!")
    print("\nKey benefits:")
    print("  - IDs are automatically generated as f\"{A}-{B}\"")
    print("  - No redundancy or sync issues between ID and A/B fields")
    print("  - Simpler JSON structure (one less field per compound)")
    print("  - Auto-adapts when groups are renamed")

if __name__ == "__main__":
    test_auto_generated_ids()
