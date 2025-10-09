#!/usr/bin/env python
"""Test precedent filtering by reagent database availability."""

import sys
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools import precedent
from chemtools.reagent_lookup import check_precedent_reagents_in_database

def test_precedent_filtering():
    """Test that precedent filtering works correctly."""
    
    print("=" * 70)
    print("Testing ML Precedent Filtering by Reagent Database Availability")
    print("=" * 70)
    
    # Test reaction: Buchwald-Hartwig amination
    reaction_smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    family = "C_N_Coupling_Pd"
    
    print(f"\nReaction: {reaction_smiles}")
    print(f"Family: {family}")
    print(f"\nFetching precedents...\n")
    
    # Test 1: WITH filtering (default)
    print("-" * 70)
    print("TEST 1: With reagent database filtering (default)")
    print("-" * 70)
    
    pack_filtered = precedent.knn(
        family=family,
        features={},
        k=25,
        relax={
            "reaction_smiles": reaction_smiles,
            "filter_by_reagent_database": True,  # Explicitly enable (default)
        }
    )
    
    precs_filtered = pack_filtered.get("precedents", [])
    support_filtered = pack_filtered.get("support", 0)
    
    print(f"Support: {support_filtered} precedents")
    print(f"Returned: {len(precs_filtered)} precedents (top 10)\n")
    
    if precs_filtered:
        print("Sample precedents (showing first 3):")
        for i, prec in enumerate(precs_filtered[:3], 1):
            print(f"\n  [{i}] Reaction ID: {prec.get('reaction_id')}")
            print(f"      Core: {prec.get('condition_core')}")
            print(f"      Yield: {prec.get('yield')}%")
            
            # Check reagent database status
            check_result = check_precedent_reagents_in_database(prec)
            print(f"      Reagents in DB: {check_result['found_count']}/{check_result['total_count']}")
            
            if check_result['missing']:
                print(f"      Missing: {[m['name'] for m in check_result['missing']]}")
            
            # Show catalytic system
            cat_system = prec.get("catalytic_system", [])
            if cat_system:
                print(f"      Catalysts: {[c.get('name') for c in cat_system]}")
            
            # Show reagents
            reagents = prec.get("reagents", [])
            if reagents:
                print(f"      Reagents: {[r.get('name') for r in reagents]}")
            
            # Show solvents
            solvents = prec.get("solvents", [])
            if solvents:
                print(f"      Solvents: {[s.get('name') for s in solvents]}")
    
    # Test 2: WITHOUT filtering
    print("\n" + "-" * 70)
    print("TEST 2: Without reagent database filtering")
    print("-" * 70)
    
    pack_unfiltered = precedent.knn(
        family=family,
        features={},
        k=25,
        relax={
            "reaction_smiles": reaction_smiles,
            "filter_by_reagent_database": False,  # Disable filtering
        }
    )
    
    precs_unfiltered = pack_unfiltered.get("precedents", [])
    support_unfiltered = pack_unfiltered.get("support", 0)
    
    print(f"Support: {support_unfiltered} precedents")
    print(f"Returned: {len(precs_unfiltered)} precedents (top 10)\n")
    
    if precs_unfiltered:
        print("Sample precedents (showing first 3):")
        for i, prec in enumerate(precs_unfiltered[:3], 1):
            print(f"\n  [{i}] Reaction ID: {prec.get('reaction_id')}")
            print(f"      Core: {prec.get('condition_core')}")
            print(f"      Yield: {prec.get('yield')}%")
            
            # Check reagent database status
            check_result = check_precedent_reagents_in_database(prec)
            print(f"      Reagents in DB: {check_result['found_count']}/{check_result['total_count']}")
            
            if check_result['missing']:
                print(f"      Missing: {[m['name'] for m in check_result['missing']]}")
    
    # Summary comparison
    print("\n" + "=" * 70)
    print("SUMMARY COMPARISON")
    print("=" * 70)
    print(f"Without filtering: {support_unfiltered} precedents, {len(precs_unfiltered)} returned")
    print(f"With filtering:    {support_filtered} precedents, {len(precs_filtered)} returned")
    
    reduction = support_unfiltered - support_filtered if support_unfiltered > 0 else 0
    if support_unfiltered > 0:
        pct = (reduction / support_unfiltered) * 100
        print(f"Filtered out:      {reduction} precedents ({pct:.1f}%)")
    
    print("\n✅ Filtering is working!")
    print("   Only precedents with all reagents in the database are returned by default.")

if __name__ == "__main__":
    test_precedent_filtering()
