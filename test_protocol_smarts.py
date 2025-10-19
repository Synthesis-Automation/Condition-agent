"""
Test script for protocol SMARTS-based recommendation

This script tests the new SMARTS matching functionality in the protocol
recommendation system.
"""

from chemtools.protocol.recommend import ProtocolRecommender, match_reaction_smarts

# Test the SMARTS matching function directly
def test_smarts_matching():
    """Test the reaction SMARTS matching function"""
    print("=" * 80)
    print("Testing SMARTS Matching Function")
    print("=" * 80)
    
    # Test case 1: Simple Suzuki coupling
    reaction = "Brc1ccc(F)cc1.OB(O)c1ccccc1>>Fc1ccc(-c2ccccc2)cc1"
    smarts = ["Br[a].OB(O)[a]>>[a][a]"]
    
    print(f"\nTest 1: Suzuki Coupling")
    print(f"Reaction: {reaction}")
    print(f"SMARTS: {smarts}")
    result = match_reaction_smarts(reaction, smarts)
    print(f"Match: {result}")
    
    # Test case 2: Alpha-arylation (from the protocol we saw)
    reaction2 = "O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1"
    smarts2 = ["O=C1CCCC1.Br[a]>>[a]C2C(CCC2)=O"]
    
    print(f"\nTest 2: Alpha-Arylation")
    print(f"Reaction: {reaction2}")
    print(f"SMARTS: {smarts2}")
    result2 = match_reaction_smarts(reaction2, smarts2)
    print(f"Match: {result2}")
    
    # Test case 3: Non-matching reaction
    reaction3 = "CCO.CC(=O)O>>CC(=O)OCC"  # Esterification
    smarts3 = ["Br[a].OB(O)[a]>>[a][a]"]  # Suzuki SMARTS
    
    print(f"\nTest 3: Non-Matching (Esterification vs Suzuki SMARTS)")
    print(f"Reaction: {reaction3}")
    print(f"SMARTS: {smarts3}")
    result3 = match_reaction_smarts(reaction3, smarts3)
    print(f"Match: {result3}")


def test_protocol_recommendation():
    """Test the full protocol recommendation with SMARTS filtering"""
    print("\n" + "=" * 80)
    print("Testing Protocol Recommendation with SMARTS Filtering")
    print("=" * 80)
    
    try:
        # Initialize recommender
        print("\nInitializing ProtocolRecommender...")
        recommender = ProtocolRecommender()
        print(f"Loaded {len(recommender.indexer.records)} protocols")
        
        # Test reaction: Alpha-arylation (matches the alpha_arylation_Pd_enamine_Dong protocol)
        reaction = "O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1"
        
        print(f"\nQuery reaction (Alpha-arylation): {reaction}")
        
        # Test with SMARTS filtering enabled
        print("\n--- With SMARTS filtering (enabled) ---")
        results_with = recommender.recommend(
            reaction_smiles=reaction,
            k=3,
            use_smarts_filter=True,
            use_standard_format=True
        )
        
        print(f"Status: {results_with['meta']['status']}")
        print(f"Num matches: {len(results_with.get('recommended_conditions', []))}")
        
        if results_with.get('recommended_conditions'):
            for i, rec in enumerate(results_with['recommended_conditions'][:3], 1):
                meta = rec.get('protocol_metadata', {})
                print(f"\n{i}. {meta.get('title', 'N/A')}")
                print(f"   Family: {meta.get('reaction_family', 'N/A')}")
                print(f"   Confidence: {rec['confidence']:.4f}")
                print(f"   SMARTS: {meta.get('reaction_smarts', [])}")
                print(f"   Chemicals: {len(rec.get('chemicals', []))} reagents")
                if rec.get('conditions'):
                    print(f"   Temperature: {rec['conditions'].get('temperature_C')}°C")
                    print(f"   Time: {rec['conditions'].get('time_h')}h")
        
        # Test with SMARTS filtering disabled
        print("\n--- With SMARTS filtering (disabled) ---")
        results_without = recommender.recommend(
            reaction_smiles=reaction,
            k=3,
            use_smarts_filter=False,
            use_standard_format=True
        )
        
        print(f"Status: {results_without['meta']['status']}")
        print(f"Num matches: {len(results_without.get('recommended_conditions', []))}")
        
        if results_without.get('recommended_conditions'):
            for i, rec in enumerate(results_without['recommended_conditions'][:3], 1):
                meta = rec.get('protocol_metadata', {})
                print(f"\n{i}. {meta.get('title', 'N/A')}")
                print(f"   Family: {meta.get('reaction_family', 'N/A')}")
                print(f"   Confidence: {rec['confidence']:.4f}")
        
        # Compare results
        print("\n--- Comparison ---")
        print(f"With SMARTS: {len(results_with.get('recommended_conditions', []))} recommendations")
        print(f"Without SMARTS: {len(results_without.get('recommended_conditions', []))} recommendations")
        
    except Exception as e:
        print(f"\nError during protocol recommendation test: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Run tests
    test_smarts_matching()
    print("\n")
    test_protocol_recommendation()
    
    print("\n" + "=" * 80)
    print("Testing Complete!")
    print("=" * 80)
