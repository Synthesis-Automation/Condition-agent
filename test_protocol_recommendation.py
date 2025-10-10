"""
Test Protocol Recommendation System

This script tests the DRFP-based protocol recommendation functionality.

Usage:
    python test_protocol_recommendation.py
"""

from chemtools.protocol import ProtocolRecommender

def test_suzuki_recommendation():
    """Test finding Suzuki protocols"""
    print("=" * 70)
    print("Test 1: Suzuki Coupling Recommendation")
    print("=" * 70)
    print()
    
    # Initialize recommender
    recommender = ProtocolRecommender()
    
    # Test query: Simple Suzuki coupling
    query = 'BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1'
    
    print(f"Query reaction: {query}")
    print()
    
    # Get recommendations
    results = recommender.recommend(
        reaction_smiles=query,
        k=5
    )
    
    print(f"Found {len(results['matches'])} matches")
    print(f"Searched {results['metadata']['num_candidates']} candidates")
    print()
    
    # Show top matches
    for i, match in enumerate(results['matches'], 1):
        print(f"{i}. Similarity: {match['similarity']:.4f}")
        print(f"   File: {match['filename']}")
        print(f"   Title: {match['source_title']}")
        print(f"   Family: {match['reaction_family']}")
        print(f"   DOI: {match['source_doi']}")
        print()


def test_borylation_recommendation():
    """Test finding borylation protocols"""
    print("=" * 70)
    print("Test 2: Borylation Recommendation")
    print("=" * 70)
    print()
    
    recommender = ProtocolRecommender()
    
    # Test query: Alkyl iodide borylation
    query = 'CCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCB1OC(C)(C)C(C)(C)O1'
    
    print(f"Query reaction: {query}")
    print()
    
    results = recommender.recommend(
        reaction_smiles=query,
        k=3
    )
    
    print(f"Top {len(results['matches'])} matches:")
    print()
    
    for i, match in enumerate(results['matches'], 1):
        print(f"{i}. Similarity: {match['similarity']:.4f}")
        print(f"   {match['source_title']}")
        print(f"   Tags: {', '.join(match['tags'][:5])}")
        print()


def test_filtered_recommendation():
    """Test recommendation with family filter"""
    print("=" * 70)
    print("Test 3: Filtered Recommendation (Palladium reactions)")
    print("=" * 70)
    print()
    
    recommender = ProtocolRecommender()
    
    # Generic aryl bromide coupling
    query = 'Brc1ccccc1.c1ccccc1>>c1ccccc1c1ccccc1'
    
    print(f"Query reaction: {query}")
    print(f"Filter: tags=['palladium']")
    print()
    
    results = recommender.recommend(
        reaction_smiles=query,
        k=5,
        tags=['palladium']
    )
    
    print(f"Found {len(results['matches'])} palladium-catalyzed protocols:")
    print()
    
    for i, match in enumerate(results['matches'], 1):
        print(f"{i}. Similarity: {match['similarity']:.4f}")
        print(f"   {match['reaction_family']}")
        print(f"   {match['source_title'][:60]}...")
        print()


def test_with_conditions():
    """Test recommendation with condition extraction"""
    print("=" * 70)
    print("Test 4: Recommendation with Condition Extraction")
    print("=" * 70)
    print()
    
    recommender = ProtocolRecommender()
    
    # Suzuki coupling
    query = 'CCBr.c1ccccc1B(O)O>>CCc1ccccc1'
    
    print(f"Query reaction: {query}")
    print()
    
    results = recommender.recommend_with_details(
        reaction_smiles=query,
        k=3
    )
    
    print(f"Top 3 matches with extracted conditions:")
    print()
    
    for i, match in enumerate(results['matches'], 1):
        print(f"{i}. {match['source_title'][:50]}...")
        print(f"   Similarity: {match['similarity']:.4f}")
        
        if 'conditions' in match:
            cond = match['conditions']
            print(f"   Catalyst: {cond.get('catalyst', 'N/A')}")
            print(f"   Ligand: {cond.get('ligand', 'N/A')}")
            print(f"   Base: {cond.get('base', 'N/A')}")
            print(f"   Solvent: {cond.get('solvent', 'N/A')}")
            print(f"   Temperature: {cond.get('temperature_C', 'N/A')} °C")
            print(f"   Time: {cond.get('time_h', 'N/A')} h")
        print()


def main():
    """Run all tests"""
    print()
    print("🧪 Testing Protocol Recommendation System")
    print()
    
    try:
        test_suzuki_recommendation()
        test_borylation_recommendation()
        test_filtered_recommendation()
        test_with_conditions()
        
        print("=" * 70)
        print("✅ All tests completed successfully!")
        print("=" * 70)
        
    except Exception as e:
        print()
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
