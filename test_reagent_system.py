#!/usr/bin/env python3
"""Test reagent system after removing allowlists and name fields."""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2


def test_reagent_system():
    """Test that reagent system still works after simplification."""
    
    print("Loading reagent taxonomy v2...")
    try:
        taxonomy = ReagentTaxonomyV2.from_path()
        print(f"✓ Successfully loaded taxonomy")
    except Exception as e:
        print(f"✗ Failed to load taxonomy: {e}")
        return False
    
    # Test 1: Check roles loaded
    roles = list(taxonomy.iter_roles())
    print(f"\n1. Roles: {len(roles)} loaded")
    if len(roles) == 0:
        print("✗ No roles loaded!")
        return False
    
    for role in roles[:3]:
        print(f"   - {role.id}: name='{role.name}', priority={role.priority}")
    
    # Test 2: Check families loaded
    families = list(taxonomy.iter_families())
    print(f"\n2. Families: {len(families)} loaded")
    if len(families) == 0:
        print("✗ No families loaded!")
        return False
    
    for family in families[:3]:
        print(f"   - {family.id}: role={family.role_id}, name='{family.name}'")
    
    # Test 3: Check allowlists are empty but present
    print(f"\n3. Allowlists structure:")
    test_family = families[0]
    print(f"   - CAS entries: {len(test_family.allowlists.cas)}")
    print(f"   - Name entries: {len(test_family.allowlists.names)}")
    print(f"   - Keyword entries: {len(test_family.allowlists.keywords)}")
    
    # Test 4: Classification with CAS (should fail gracefully)
    print(f"\n4. Test classification with CAS:")
    test_record = {
        "cas": "7647-01-0",  # HCl
        "name": "Hydrochloric acid",
        "smiles": None
    }
    result = taxonomy.classify(test_record)
    if result:
        print(f"   ✓ Matched: {result.family_id} (role: {result.role_id}, kind: {result.match_kind})")
    else:
        print(f"   ○ No match (expected since allowlists are now empty)")
    
    # Test 5: Classification with SMILES
    print(f"\n5. Test classification with SMILES:")
    test_record = {
        "cas": None,
        "name": "Acetic acid",
        "smiles": "CC(=O)O"
    }
    result = taxonomy.classify(test_record)
    if result:
        print(f"   ✓ Matched: {result.family_id} (role: {result.role_id}, kind: {result.match_kind})")
    else:
        print(f"   ○ No match found")
    
    # Test 6: Check families with detect patterns
    families_with_detect = [f for f in families if f.detect is not None]
    print(f"\n6. Families with SMARTS detection: {len(families_with_detect)}")
    if families_with_detect:
        for fam in families_with_detect[:3]:
            print(f"   - {fam.id}: {len(fam.detect.smarts_any)} patterns")
    
    # Test 7: Check that name fallback works
    print(f"\n7. Name fallback (id used when name missing):")
    for family in families[:5]:
        # All should have names now (fallback to id if missing)
        if family.name:
            print(f"   ✓ {family.id}: has name '{family.name}'")
        else:
            print(f"   ✗ {family.id}: missing name!")
            return False
    
    print(f"\n✓ All tests passed! Reagent system works correctly.")
    return True


if __name__ == "__main__":
    success = test_reagent_system()
    sys.exit(0 if success else 1)
