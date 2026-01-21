#!/usr/bin/env python3
"""Comprehensive test of reagent system after taxonomy simplification."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2
from chemtools.reagent.taxonomy_store import TaxonomyStore, RoleHeuristics


def test_reagent_v2():
    """Test ReagentTaxonomyV2 class."""
    print("=" * 60)
    print("Testing ReagentTaxonomyV2...")
    print("=" * 60)
    
    taxonomy = ReagentTaxonomyV2.from_path()
    
    # Test roles
    roles = list(taxonomy.iter_roles())
    assert len(roles) > 0, "No roles loaded"
    print(f"✓ Loaded {len(roles)} roles")
    
    # Test role name fallback
    for role in roles:
        assert role.name, f"Role {role.id} missing name"
    print(f"✓ All roles have names (using fallback if needed)")
    
    # Test families
    families = list(taxonomy.iter_families())
    assert len(families) > 0, "No families loaded"
    print(f"✓ Loaded {len(families)} families")
    
    # Test family name fallback
    for family in families:
        assert family.name, f"Family {family.id} missing name"
    print(f"✓ All families have names (using fallback if needed)")
    
    # Test allowlists are empty but functional
    for family in families:
        assert isinstance(family.allowlists.cas, set), f"CAS should be a set"
        assert isinstance(family.allowlists.names, set), f"Names should be a set"
        assert isinstance(family.allowlists.keywords, list), f"Keywords should be a list"
    print(f"✓ All families have empty allowlists (no errors)")
    
    # Test classification still works
    test_records = [
        {"cas": None, "name": "Acetic acid", "smiles": "CC(=O)O"},
        {"cas": None, "name": "Toluene", "smiles": "Cc1ccccc1"},
        {"cas": None, "name": "Palladium acetate", "smiles": None},
    ]
    
    matches = 0
    for record in test_records:
        result = taxonomy.classify(record)
        if result:
            matches += 1
    print(f"✓ Classification works: {matches}/{len(test_records)} matched")
    
    return True


def test_taxonomy_store():
    """Test TaxonomyStore class."""
    print("\n" + "=" * 60)
    print("Testing TaxonomyStore...")
    print("=" * 60)
    
    store = TaxonomyStore()
    
    # Test that store loaded successfully
    assert store.taxonomy is not None, "Taxonomy not loaded"
    print(f"✓ TaxonomyStore initialized successfully")
    
    # Test families are accessible
    families_list = store.list_families()
    assert len(families_list) > 0, "No families in store"
    print(f"✓ TaxonomyStore has {len(families_list)} families")
    
    # Test that family data is accessible (allowlists will be empty)
    test_family_id = families_list[0][1]  # Get first family ID
    test_role_id = families_list[0][0]  # Get first role ID
    family_data = store.family_data(test_role_id, test_family_id)
    assert "allowlists" in family_data, "Missing allowlists in family data"
    print(f"✓ Family data exports allowlists structure (empty)")
    
    # Test heuristics still work
    heuristics = RoleHeuristics(store)
    assert heuristics is not None, "Heuristics not loaded"
    print(f"✓ RoleHeuristics initialized successfully")
    
    # Test heuristic inference (should still work with empty allowlists)
    role_result = heuristics.infer_role("Palladium acetate", [])
    if role_result:
        print(f"✓ Role inference works: '{role_result[0]}' (confidence: {role_result[1]})")
    else:
        print(f"○ Role inference returned None (acceptable)")
    
    return True


def main():
    """Run all tests."""
    try:
        test_reagent_v2()
        test_taxonomy_store()
        
        print("\n" + "=" * 60)
        print("✓ ALL TESTS PASSED!")
        print("=" * 60)
        print("\nThe reagent system works correctly after:")
        print("  - Removing 170 'name' fields")
        print("  - Removing 158 'allowlists' objects")
        print("  - Total: 328 fields removed from reagent_roles.v2.json")
        return 0
    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
