"""
Demo script for UnifiedRecommender - DRFP-based unified recommendation system.

This script demonstrates the unified recommendation engine that combines
protocols and rules into a single similarity-based recommendation system.
"""

from chemtools.recommend.unified import UnifiedRecommender


def demo_unified_recommender():
    """Demonstrate unified recommendation with various query reactions."""
    
    print("="*80)
    print("UnifiedRecommender Demo")
    print("="*80)
    
    # Initialize recommender
    print("\n[1] Loading unified index...")
    recommender = UnifiedRecommender("build/unified_index_complete")
    
    # Show statistics
    stats = recommender.get_statistics()
    print(f"\nIndex Statistics:")
    print(f"  Version: {stats.get('version', stats.get('build_info', {}).get('version', 'unknown'))}")
    print(f"  Protocols: {stats.get('num_protocols', stats.get('protocols', {}).get('count', 0))}")
    print(f"  Rules: {stats.get('num_rules', stats.get('rules', {}).get('count', 0))}")
    print(f"  Fingerprints: {stats.get('num_fingerprints', stats.get('drfp', {}).get('computed', 0))}")
    print(f"  Families: {stats.get('num_families', 'N/A')}")
    
    # Test queries (more realistic reactions with functional groups)
    test_reactions = [
        {
            "name": "Suzuki coupling (Br-aryl + boronic acid)",
            "smiles": "CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1",
            "description": "Suzuki coupling of 4-bromoacetophenone with phenylboronic acid"
        },
        {
            "name": "Buchwald-Hartwig amination (Br-aryl + amine)",
            "smiles": "Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1",
            "description": "C-N coupling of bromobenzene with benzylamine"
        },
        {
            "name": "RCM (diene cyclization)",
            "smiles": "C=CCN(C(=O)OC(C)(C)C)CC=C>>CC(C)(C)OC(=O)N1CCC=C1",
            "description": "Ring-closing metathesis of N-Boc-diallylamine"
        },
        {
            "name": "Sonogashira coupling (I-aryl + terminal alkyne)",
            "smiles": "Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
            "description": "Pd/Cu-catalyzed coupling of iodobenzene with phenylacetylene"
        },
        {
            "name": "Amide formation (acid + amine)",
            "smiles": "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1",
            "description": "Coupling of benzoic acid with benzylamine"
        },
    ]
    
    # Run recommendations for each query
    for i, query in enumerate(test_reactions, start=1):
        print(f"\n{'='*80}")
        print(f"[{i+1}] Query: {query['name']}")
        print(f"Description: {query['description']}")
        print(f"SMILES: {query['smiles']}")
        print(f"{'='*80}")
        
        try:
            # Get top 5 recommendations
            results = recommender.recommend(
                query['smiles'],
                top_k=5,
                min_similarity=0.3
            )
            
            if not results:
                print("  No recommendations found (similarity < 0.3)")
                continue
            
            print(f"\n  Top {len(results)} Recommendations:")
            print(f"  {'Rank':<6} {'Type':<10} {'Similarity':<12} {'Name':<50}")
            print(f"  {'-'*78}")
            
            for result in results:
                source_type_icon = "📋" if result.source_type == "protocol" else "📖"
                print(f"  {result.rank:<6} {source_type_icon} {result.source_type:<8} "
                      f"{result.similarity:>6.3f}      {result.name[:48]}")
            
            # Show details for top match
            top_result = results[0]
            print(f"\n  🏆 Best Match Details:")
            print(f"     ID: {top_result.id}")
            print(f"     Family: {top_result.family}")
            print(f"     Version: {top_result.version}")
            print(f"     Tags: {', '.join(top_result.tags[:5])}")
            
            # Test loading full details
            print(f"\n  Loading full source data...")
            full_data = recommender.get_source_details(top_result.id)
            if full_data:
                print(f"     ✅ Successfully loaded {len(str(full_data))} characters of data")
            else:
                print(f"     ⚠️  Could not load full data")
        
        except Exception as e:
            print(f"  ❌ Error: {e}")
            import traceback
            traceback.print_exc()
    
    # Test filtering by source type
    print(f"\n{'='*80}")
    print("[Additional Tests] Filtering by Source Type")
    print(f"{'='*80}")
    
    test_smiles = "c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    
    # Protocols only
    print("\n  Protocols only (top 3):")
    protocol_results = recommender.recommend(
        test_smiles,
        top_k=3,
        source_types=["protocol"]
    )
    for r in protocol_results:
        print(f"    • {r.name[:50]} (sim: {r.similarity:.3f})")
    
    # Rules only
    print("\n  Rules only (top 3):")
    rule_results = recommender.recommend(
        test_smiles,
        top_k=3,
        source_types=["rule"]
    )
    for r in rule_results:
        print(f"    • {r.name[:50]} (sim: {r.similarity:.3f})")
    
    print(f"\n{'='*80}")
    print("Demo Complete!")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    demo_unified_recommender()
