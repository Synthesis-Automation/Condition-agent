"""
Quick debug script to test fusion recommendation directly.
"""

from chemtools.ml.fusion_recommender import recommend_with_fusion

reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
family = "C_N_Coupling_Pd"

print("Running fusion recommendation...")
try:
    results = recommend_with_fusion(reaction, family, k=50, top_n=5)
    print("✅ Fusion recommendation successful!")
    print(f"\nReturned keys: {results.keys()}")
    
    # Check evidence structure
    evidence = results.get('evidence', {})
    print(f"\nEvidence keys: {evidence.keys() if evidence else 'None'}")
    
    if 'precedents' in evidence:
        prec = evidence['precedents']
        print(f"\nPrecedent evidence keys: {prec.keys()}")
        print(f"  Coverage: {prec.get('coverage')}")
        print(f"  Diversity: {prec.get('diversity_score')}")
        print(f"  Num reactions: {len(prec.get('reactions', []))}")
    
    if 'analytics' in evidence:
        analytics = evidence['analytics']
        print(f"\nAnalytics keys: {analytics.keys()}")
        print(f"  Dataset size: {analytics.get('dataset_size')}")
    
    # Check adaptive_weights
    aw = results.get('adaptive_weights', {})
    print(f"\nAdaptive weights keys: {aw.keys()}")
    print(f"  Weights: {aw.get('weights')}")
    print(f"  Reasoning: {aw.get('reasoning')}")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
