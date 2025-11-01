#!/usr/bin/env python3
"""
Simple debug test to understand cross-family behavior.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))


def debug_cross_family():
    """Debug the cross-family functionality step by step."""
    
    print("🔍 Debug Cross-Family Functionality")
    print("=" * 50)
    
    test_reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    print(f"Test Reaction: {test_reaction}")
    print()
    
    try:
        from chemtools import chem
        
        # Test family detection first
        print("1️⃣ Testing Family Detection:")
        from chemtools.detection import detect_reaction
        detection_result = detect_reaction(test_reaction)
        print(f"Detection result: {detection_result}")
        print()
        
        # Test direct API
        print("2️⃣ Testing Direct Recommendation API:")
        result = chem.recommend.from_reaction(
            reaction=test_reaction,
            k=10,
            relax={
                'search_all_families': True,
                'reaction_type_threshold': 0.10,
                'mechanism_similarity_threshold': 0.30,
                'mechanism_weight': 0.25,
            }
        )
        
        print(f"Raw result keys: {list(result.keys())}")
        recommendations = result.get('recommendations', [])
        print(f"Found {len(recommendations)} recommendations")
        
        if recommendations:
            first_rec = recommendations[0]
            print(f"First recommendation keys: {list(first_rec.keys())}")
            
            conditions = first_rec.get('conditions', {})
            print(f"Conditions keys: {list(conditions.keys())}")
            
            cross_meta = conditions.get('cross_family_metadata', {})
            print(f"Cross-family metadata: {cross_meta}")
        
        print()
        
        # Test structured API
        print("3️⃣ Testing Structured API:")
        structured_result = chem.recommend.conditions(
            reaction=test_reaction,
            k=10,
            search_all_families=True,
            relax={
                'reaction_type_threshold': 0.10,
                'mechanism_similarity_threshold': 0.30,
                'mechanism_weight': 0.25,
            }
        )
        
        print(f"Structured result keys: {list(structured_result.keys())}")
        
        detection = structured_result.get('detection', {})
        print(f"Detection info: {detection}")
        
        recs = structured_result.get('recommendations', [])
        print(f"Structured recommendations count: {len(recs)}")
        
        if recs:
            print("\n4️⃣ Analyzing First Recommendation:")
            first_rec = recs[0]
            print(f"Recommendation keys: {list(first_rec.keys())}")
            
            conditions = first_rec.get('conditions', {})
            print(f"Conditions type: {type(conditions)}")
            if isinstance(conditions, dict):
                print(f"Conditions keys: {list(conditions.keys())}")
            cross_meta = conditions.get('cross_family_metadata')
            print(f"Cross-family metadata: {cross_meta}")
            
            # Check the raw precedents
            print("\n5️⃣ Checking Raw Precedents:")
            precedents = structured_result.get('precedents', [])
            print(f"Precedents type: {type(precedents)}")
            print(f"Total precedents: {len(precedents) if isinstance(precedents, (list, dict)) else 'N/A'}")
            
            # Also check precedents_used
            precedents_used = structured_result.get('precedents_used', [])
            print(f"Precedents used: {len(precedents_used) if isinstance(precedents_used, list) else 'N/A'}")
            
            if precedents_used and isinstance(precedents_used, list):
                first_prec = precedents_used[0]
                print(f"First precedent_used keys: {list(first_prec.keys())}")
                print(f"Has mechanism_similarity: {'mechanism_similarity' in first_prec}")
                print(f"Has cross_family_metadata: {'cross_family_metadata' in first_prec}")
                if 'cross_family_metadata' in first_prec:
                    print(f"Metadata content: {first_prec['cross_family_metadata']}")
                if 'rxn_type' in first_prec:
                    print(f"Reaction type: {first_prec['rxn_type']}")
        
    except Exception as e:
        print(f"❌ Debug failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    debug_cross_family()