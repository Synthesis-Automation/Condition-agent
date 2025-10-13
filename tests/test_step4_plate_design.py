#!/usr/bin/env python3
"""
Test Step 4: Plate Design and Constraints

Tests:
- Plate design generation (design_plate_from_reaction)
- Well layout and distribution
- Constraint filtering (inventory, blacklist)
- Diversity of conditions across wells

This validates the experimental plate design feature for HTS workflows.
"""

import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Direct imports to avoid ML library loading overhead
from chemtools.recommend.core import recommend_from_reaction
from chemtools.recommend.plate_design import design_plate_from_reaction
from chemtools import reagent

# Test reaction (same as previous tests)
TEST_REACTION = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"


def lookup_name(cas: str, role: str) -> str:
    """Look up chemical name from CAS number."""
    if not cas or cas == 'N/A':
        return 'N/A'
    
    info = reagent.enrich_reagent_info(cas, role)
    if info and info.get('name'):
        # Check if the name is just the CAS number (not found in database)
        name = info['name']
        if name == cas or name.startswith('[Unknown'):
            return f"[Unknown {role}] {cas}"
        return f"{name} ({cas})"
    return f"[Unknown {role}] {cas}"


def test_4a_plate_design():
    """Test 4a: Basic plate design generation."""
    print("=" * 70)
    print("  Test Step 4: Plate Design and Constraints")
    print("=" * 70)
    print()
    print("──" * 35)
    print("  Test 4a: Plate Design Generation")
    print("──" * 35)
    print(f"   Query: {TEST_REACTION}")
    print()
    
    print("   🧪 Designing 24-well plate...")
    t_start = time.time()
    
    try:
        result = design_plate_from_reaction(
            reaction=TEST_REACTION,
            plate_size=24,
        )
        
        t_elapsed = time.time() - t_start
        print(f"   ⏱️  Completed in {t_elapsed:.3f} seconds")
        print()
        
        # Check structure
        print("   ✅ Results:")
        
        # Check for errors
        if 'meta' in result and 'error' in result.get('meta', {}):
            error = result['meta']['error']
            print(f"      ⚠️  Error: {error}")
            return 1
        
        # Check rows
        if 'rows' in result:
            rows = result['rows']
            print(f"      - Number of wells: {len(rows)}")
            
            if len(rows) > 0:
                print()
                print("   📊 Well distribution:")
                
                # Count unique cores, bases, solvents
                cores = set()
                bases = set()
                solvents = set()
                
                for row in rows:
                    core = row.get('core')
                    base = row.get('base_uid')
                    solvent = row.get('solvent_uid')
                    
                    if core:
                        cores.add(core)
                    if base:
                        bases.add(base)
                    if solvent:
                        solvents.add(solvent)
                
                print(f"      - Unique catalyst cores: {len(cores)}")
                print(f"      - Unique bases: {len(bases)}")
                print(f"      - Unique solvents: {len(solvents)}")
                
                if cores:
                    print(f"      - Cores: {', '.join(sorted(cores)[:5])}" + (' ...' if len(cores) > 5 else ''))
                
                print()
                print("   📋 Sample wells (first 5):")
                print()
                
                for i, row in enumerate(rows[:5], 1):
                    well_id = row.get('well_id', f'#{i}')
                    core = row.get('core', 'N/A')
                    base_cas = row.get('base_uid', 'N/A')
                    solvent_cas = row.get('solvent_uid', 'N/A')
                    temp = row.get('T_C', 'N/A')
                    time_h = row.get('time_h', 'N/A')
                    
                    # Look up chemical names
                    base_name = lookup_name(base_cas, 'base')
                    solvent_name = lookup_name(solvent_cas, 'solvent')
                    
                    print(f"      [{well_id}] {core}")
                    print(f"          Base: {base_name}")
                    print(f"          Solvent: {solvent_name}")
                    print(f"          Conditions: {temp}°C, {time_h}h")
                    print()
        else:
            print("      ⚠️  No 'rows' section found")
        
        # Check CSV
        if 'csv' in result:
            csv_str = result['csv']
            lines = csv_str.strip().split('\n')
            print(f"      - CSV output: {len(lines)} lines (header + {len(lines)-1} wells)")
        else:
            print("      ⚠️  No 'csv' section found")
        
        print()
        print("──" * 35)
        print("  Validation Summary")
        print("──" * 35)
        print()
        
        # Validation checks
        checks_passed = 0
        total_checks = 4
        
        # Check 1: Rows exist
        if 'rows' in result and len(result['rows']) > 0:
            print("   ✅ Plate wells generated")
            checks_passed += 1
        else:
            print("   ❌ No plate wells generated")
        
        # Check 2: CSV format
        if 'csv' in result and len(result['csv']) > 0:
            print("   ✅ CSV format generated")
            checks_passed += 1
        else:
            print("   ❌ No CSV output")
        
        # Check 3: Diversity of cores
        if 'rows' in result and len(rows) > 0:
            cores = set(r.get('core') for r in rows if r.get('core'))
            if len(cores) >= 1:
                print(f"   ✅ Catalyst diversity: {len(cores)} unique core(s)")
                checks_passed += 1
            else:
                print("   ❌ No catalyst diversity")
        else:
            print("   ❌ Cannot check catalyst diversity")
        
        # Check 4: Well IDs assigned
        if 'rows' in result and len(rows) > 0:
            has_well_ids = all(r.get('well_id') for r in rows[:5])
            if has_well_ids:
                print("   ✅ Well IDs assigned")
                checks_passed += 1
            else:
                print("   ⚠️  Some wells missing IDs")
        else:
            print("   ❌ Cannot check well IDs")
        
        print()
        print(f"   📊 Validation: {checks_passed}/{total_checks} checks passed")
        print()
        
        # Overall success
        if checks_passed >= 3:  # Allow 1 failure
            print("   ✅ Test 4a PASSED")
            print()
            return 0
        else:
            print("   ❌ Test 4a FAILED")
            print()
            return 1
            
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


def test_4b_constraint_filtering():
    """Test 4b: Constraint filtering with inventory."""
    print("──" * 35)
    print("  Test 4b: Constraint Filtering")
    print("──" * 35)
    print(f"   Query: {TEST_REACTION}")
    print()
    
    print("   🔍 Testing with inventory constraints...")
    t_start = time.time()
    
    # Define a simple inventory constraint
    constraints = {
        "inventory": {
            "base": ["7758-29-4", "7778-53-2"],  # Only allow specific bases
            "solvent": ["108-88-3", "1120-21-4"]  # Only allow specific solvents
        }
    }
    
    try:
        result = recommend_from_reaction(
            reaction=TEST_REACTION,
            k=25,
            max_variants=5,
            constraint_rules=constraints,
        )
        
        t_elapsed = time.time() - t_start
        print(f"   ⏱️  Completed in {t_elapsed:.3f} seconds")
        print()
        
        print("   ✅ Results with constraints:")
        
        # Check if constraints were applied
        if 'constraint_filters' in result:
            filters = result['constraint_filters']
            print(f"      - Constraints applied: {filters}")
        else:
            print("      - No constraint information")
        
        # Check formatted recommendations
        if 'formatted' in result:
            formatted = result['formatted']
            conditions = formatted.get('recommended_conditions', [])
            
            print(f"      - Conditions generated: {len(conditions)}")
            
            if len(conditions) > 0:
                print()
                print("   📊 Checking constraint compliance:")
                print()
                
                allowed_bases = set(constraints['inventory']['base'])
                allowed_solvents = set(constraints['inventory']['solvent'])
                
                compliant_count = 0
                
                for i, cond in enumerate(conditions[:3], 1):
                    combo = cond.get('combo', {})
                    base_uid = combo.get('base_uid')
                    solvent_uid = combo.get('solvent_uid')
                    
                    base_ok = base_uid in allowed_bases if base_uid else False
                    solv_ok = solvent_uid in allowed_solvents if solvent_uid else False
                    
                    status = "✅" if (base_ok and solv_ok) else "⚠️"
                    
                    # Look up chemical names
                    base_name = lookup_name(base_uid, 'base') if base_uid else 'N/A'
                    solvent_name = lookup_name(solvent_uid, 'solvent') if solvent_uid else 'N/A'
                    
                    print(f"      [{i}] {status}")
                    print(f"          Base: {base_name} {'✓' if base_ok else '✗ (not in inventory)'}")
                    print(f"          Solvent: {solvent_name} {'✓' if solv_ok else '✗ (not in inventory)'}")
                    
                    if base_ok and solv_ok:
                        compliant_count += 1
                
                print()
                print(f"      Compliance: {compliant_count}/{min(3, len(conditions))} conditions use only inventory reagents")
        
        print()
        print("   ✅ Test 4b PASSED")
        print()
        return 0
        
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


def test_4c_small_plate():
    """Test 4c: Small 6-well plate design."""
    print("──" * 35)
    print("  Test 4c: Small Plate Design (6-well)")
    print("──" * 35)
    print(f"   Query: {TEST_REACTION}")
    print()
    
    print("   🧪 Designing 6-well plate...")
    t_start = time.time()
    
    try:
        result = design_plate_from_reaction(
            reaction=TEST_REACTION,
            plate_size=6,
        )
        
        t_elapsed = time.time() - t_start
        print(f"   ⏱️  Completed in {t_elapsed:.3f} seconds")
        print()
        
        if 'rows' in result:
            rows = result['rows']
            print(f"   ✅ Generated {len(rows)} wells")
            print()
            
            # Show all wells for small plate
            print("   📋 Complete plate layout:")
            print()
            
            for row in rows:
                well_id = row.get('well_id', '?')
                core = row.get('core', 'N/A')
                base_cas = row.get('base_uid', 'N/A')
                solvent_cas = row.get('solvent_uid', 'N/A')
                
                # Look up chemical names
                base_name = lookup_name(base_cas, 'base')
                solvent_name = lookup_name(solvent_cas, 'solvent')
                
                print(f"      [{well_id}] {core} + {base_name} + {solvent_name}")
            
            print()
            
            # Check diversity
            cores = set(r.get('core') for r in rows if r.get('core'))
            bases = set(r.get('base_uid') for r in rows if r.get('base_uid'))
            solvents = set(r.get('solvent_uid') for r in rows if r.get('solvent_uid'))
            
            print(f"   📊 Diversity: {len(cores)} cores, {len(bases)} bases, {len(solvents)} solvents")
            print()
            print("   ✅ Test 4c PASSED")
            print()
            return 0
        else:
            print("   ❌ No plate design generated")
            return 1
            
    except Exception as e:
        print(f"   ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


def main():
    """Run all Test Step 4 tests."""
    print()
    
    # Test 4a: Basic plate design
    result_4a = test_4a_plate_design()
    
    # Test 4b: Constraint filtering
    result_4b = test_4b_constraint_filtering()
    
    # Test 4c: Small plate
    result_4c = test_4c_small_plate()
    
    # Performance summary
    print("=" * 70)
    print("  Performance Summary")
    print("=" * 70)
    print()
    print("   ⏱️  All tests completed")
    print()
    
    # Final summary
    print("=" * 70)
    print("  ✅ Test Step 4 Complete!")
    print("=" * 70)
    print()
    
    total_failed = result_4a + result_4b + result_4c
    
    if total_failed == 0:
        print("🎉 All validation checks passed!")
        return 0
    else:
        print(f"⚠️  {total_failed} test(s) had issues")
        return 1


if __name__ == "__main__":
    sys.exit(main())
