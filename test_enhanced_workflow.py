#!/usr/bin/env python3
"""
Test script to demonstrate the enhanced cross-family workflow with marking instead of filtering.

This demonstrates the new approach where all precedents are kept but marked and ranked appropriately.
"""

import sys
import json
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_enhanced_cross_family():
    """Test the enhanced cross-family recommendation with marking."""
    
    # Test reaction: C-N coupling (should detect as C_N_Coupling_Pd)
    test_reaction = "CC(C)c1ccc(Br)cc1.Nc2ccccc2>>CC(C)c1ccc(Nc2ccccc2)cc1"
    
    print("🧪 Testing Enhanced Cross-Family Workflow")
    print(f"📝 Reaction: {test_reaction}")
    print()
    
    try:
        from chemtools import chem
        
        # Test 1: Family-specific search (baseline)
        print("1️⃣ Family-Specific Search (Baseline):")
        family_result = chem.recommend.recommend_from_reaction(
            test_reaction,
            k=20,
            search_all_families=False,
            relax={'debug_timing': True}
        )
        
        family_precs = family_result.get('recommendations', [])
        print(f"   Found {len(family_precs)} family-specific precedents")
        if family_precs:
            detected_family = family_result.get('detection', {}).get('canonical_family', 'Unknown')
            print(f"   Detected family: {detected_family}")
        print()
        
        # Test 2: Cross-family search with enhanced marking
        print("2️⃣ Enhanced Cross-Family Search (with marking):")
        cross_result = chem.recommend.recommend_from_reaction(
            test_reaction,
            k=50,
            search_all_families=True,
            relax={
                'reaction_type_threshold': 0.10,  # 10% threshold  
                'mechanism_similarity_threshold': 0.30,  # Lower threshold to keep more
                'mechanism_weight': 0.25,
                'debug_timing': True
            },
            rerank_strategy='none'  # Use mechanism reranking instead
        )
        
        cross_precs = cross_result.get('recommendations', [])
        print(f"   Found {len(cross_precs)} cross-family precedents")
        
        # Analyze the marking results
        if cross_precs:
            analyze_marking_results(cross_precs[:10])  # Analyze top 10
        print()
        
        # Test 3: Show comparison
        print("3️⃣ Comparison Summary:")
        print(f"   Family-specific: {len(family_precs)} precedents")
        print(f"   Cross-family: {len(cross_precs)} precedents")
        
        if cross_precs:
            compatible_count = sum(1 for p in cross_precs 
                                 if p.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'compatible')
            moderate_count = sum(1 for p in cross_precs 
                               if p.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'moderate')
            incompatible_count = sum(1 for p in cross_precs 
                                   if p.get('conditions', {}).get('cross_family_metadata', {}).get('mechanism_status') == 'incompatible')
            
            print(f"   Mechanism compatibility breakdown:")
            print(f"     • Compatible: {compatible_count}")
            print(f"     • Moderate: {moderate_count}")  
            print(f"     • Incompatible: {incompatible_count}")
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()


def analyze_marking_results(precedents):
    """Analyze the marking results from cross-family search."""
    
    print("   📊 Top 10 Precedents Analysis:")
    print("   " + "="*80)
    
    for i, prec in enumerate(precedents, 1):
        conditions = prec.get('conditions', {})
        cross_meta = conditions.get('cross_family_metadata', {})
        
        # Extract key info
        core = prec.get('core', 'Unknown')
        family = cross_meta.get('precedent_family', 'Unknown')
        mech_sim = cross_meta.get('mechanism_similarity', 0)
        mech_status = cross_meta.get('mechanism_status', 'unknown')
        type_status = cross_meta.get('reaction_type_status', 'unknown')
        similarity = conditions.get('similarity', 0)
        
        # Format status indicators
        mech_icon = {'compatible': '✅', 'moderate': '⚠️', 'incompatible': '❌'}.get(mech_status, '❓')
        type_icon = {'well_represented': '📊', 'underrepresented': '📉'}.get(type_status, '❓')
        
        print(f"   {i:2d}. {core[:30]:<30} | {family:<15}")
        print(f"       Sim: {similarity:.3f} | Mech: {mech_sim:.3f} {mech_icon} | Type: {type_icon}")
        print()


def test_api_compatibility():
    """Test that the API still works with the new parameters."""
    
    print("🌐 Testing API Compatibility")
    
    try:
        from chemtools.contracts import RecommendFromReactionRequest
        from app.services.recommendation_service import recommend_from_reaction
        
        # Test API request with new parameters
        request = RecommendFromReactionRequest(
            reaction="CC(C)c1ccc(Br)cc1.Nc2ccccc2>>CC(C)c1ccc(Nc2ccccc2)cc1",
            k=20,
            search_all_families=True,
            reaction_type_threshold=0.15,
            mechanism_similarity_threshold=0.4,
            mechanism_weight=0.3
        )
        
        result = recommend_from_reaction(request)
        
        print(f"   ✅ API request succeeded")
        print(f"   📊 Found {len(result.get('recommendations', []))} recommendations")
        
        # Check for cross-family metadata
        recs = result.get('recommendations', [])
        if recs:
            first_rec = recs[0]
            has_cross_meta = 'cross_family_metadata' in first_rec.get('conditions', {})
            print(f"   🏷️ Cross-family metadata present: {'✅' if has_cross_meta else '❌'}")
        
    except Exception as e:
        print(f"   ❌ API test failed: {e}")


if __name__ == "__main__":
    test_enhanced_cross_family()
    print()
    test_api_compatibility()