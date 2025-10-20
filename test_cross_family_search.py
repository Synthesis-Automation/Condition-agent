"""
Test script for cross-family ML recommendation search.

This script demonstrates the new search_all_families feature that allows
searching across all reaction datasets without filtering by reaction type.

Usage:
    python test_cross_family_search.py
"""

import sys
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT))

from chemtools import chem


def test_cross_family_search():
    """Test cross-family search vs standard family-specific search."""
    
    # Test reaction: C-N coupling (Buchwald-Hartwig style)
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print("=" * 80)
    print("Cross-Family ML Recommendation Search Test")
    print("=" * 80)
    print(f"\nReaction: {reaction}")
    print()
    print("NOTE: Cross-family search currently uses feature-based similarity")
    print("      (DRFP disabled by default for performance)")
    print("      Results may have lower confidence but are much faster.")
    print()
    
    # Test 1: Standard family-specific search (auto-detect)
    print("-" * 80)
    print("Test 1: STANDARD SEARCH (family-specific, auto-detect)")
    print("-" * 80)
    
    result_standard = chem.recommend.conditions(
        reaction=reaction,
        k=50,
        limit=3,
        search_all_families=False
    )
    
    detected_family = result_standard.get('detection', {}).get('detected_reaction_type', 'Unknown')
    recs_standard = result_standard.get('recommendations', [])
    
    print(f"Detected family: {detected_family}")
    print(f"Search mode: {result_standard.get('detection', {}).get('search_mode', 'unknown')}")
    print(f"Number of recommendations: {len(recs_standard)}")
    
    if recs_standard:
        print("\nTop recommendation:")
        rec = recs_standard[0]
        chemicals = rec.get('chemicals', [])
        conditions = rec.get('conditions', {})
        confidence = rec.get('confidence', 0)
        
        # Show catalyst/core
        catalysts = [c for c in chemicals if c.get('role') in ['catalyst', 'ligand', 'metal_precursor']]
        if catalysts:
            print(f"  Catalyst system: {', '.join(c.get('name', 'unknown') for c in catalysts)}")
        
        # Show base and solvent
        base = next((c for c in chemicals if c.get('role') == 'base'), None)
        solvent = next((c for c in chemicals if c.get('role') == 'solvent'), None)
        
        if base:
            print(f"  Base: {base.get('name', 'unknown')}")
        if solvent:
            print(f"  Solvent: {solvent.get('name', 'unknown')}")
        
        # Show conditions
        if conditions.get('temperature'):
            print(f"  Temperature: {conditions['temperature'].get('value')}°C")
        if conditions.get('time'):
            print(f"  Time: {conditions['time'].get('value')}h")
        
        print(f"  Confidence: {confidence:.2f}")
    
    print()
    
    # Test 2: Cross-family search (all datasets)
    print("-" * 80)
    print("Test 2: CROSS-FAMILY SEARCH (search all datasets)")
    print("-" * 80)
    
    result_cross = chem.recommend.conditions(
        reaction=reaction,
        k=100,  # Use more precedents for broader search
        limit=3,
        search_all_families=True
    )
    
    detected_family_cross = result_cross.get('detection', {}).get('detected_reaction_type', 'All')
    recs_cross = result_cross.get('recommendations', [])
    
    print(f"Detected family: {detected_family_cross}")
    print(f"Search mode: {result_cross.get('detection', {}).get('search_mode', 'unknown')}")
    print(f"Number of recommendations: {len(recs_cross)}")
    
    if recs_cross:
        print("\nTop recommendation:")
        rec = recs_cross[0]
        chemicals = rec.get('chemicals', [])
        conditions = rec.get('conditions', {})
        confidence = rec.get('confidence', 0)
        
        # Show catalyst/core
        catalysts = [c for c in chemicals if c.get('role') in ['catalyst', 'ligand', 'metal_precursor']]
        if catalysts:
            print(f"  Catalyst system: {', '.join(c.get('name', 'unknown') for c in catalysts)}")
        
        # Show base and solvent
        base = next((c for c in chemicals if c.get('role') == 'base'), None)
        solvent = next((c for c in chemicals if c.get('role') == 'solvent'), None)
        
        if base:
            print(f"  Base: {base.get('name', 'unknown')}")
        if solvent:
            print(f"  Solvent: {solvent.get('name', 'unknown')}")
        
        # Show conditions
        if conditions.get('temperature'):
            print(f"  Temperature: {conditions['temperature'].get('value')}°C")
        if conditions.get('time'):
            print(f"  Time: {conditions['time'].get('value')}h")
        
        print(f"  Confidence: {confidence:.2f}")
    
    print()
    
    # Compare results
    print("-" * 80)
    print("COMPARISON")
    print("-" * 80)
    print(f"Standard search found: {len(recs_standard)} recommendations from {detected_family} family")
    print(f"Cross-family search found: {len(recs_cross)} recommendations from ALL families")
    print()
    
    # Show which families were found in cross-family search
    if 'precedents_used' in result_cross:
        precs = result_cross['precedents_used'].get('top_precedents', [])
        if precs:
            families_found = set()
            for p in precs:
                rxn_type = p.get('rxn_type', 'Unknown')
                if rxn_type:
                    families_found.add(rxn_type)
            
            print(f"Families represented in cross-family precedents: {', '.join(sorted(families_found))}")
    
    print()
    print("=" * 80)
    print("✓ Test completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    test_cross_family_search()
