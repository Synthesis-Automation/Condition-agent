"""
Quick start guide for HTE recommendation system
"""
from chemtools.HTE import HTERecommender, format_result

print("="*80)
print("HTE RECOMMENDATION SYSTEM - QUICK START")
print("="*80)

print("\n🚀 Step 1: Initialize recommender")
print("   Code: recommender = HTERecommender()")
recommender = HTERecommender()

print("\n📊 Step 2: Check database statistics")
stats = recommender.get_statistics()
print(f"   Database: {stats['total_experiments']:,} experiments")
print(f"   Reaction types: {stats['reaction_types']}")
print(f"   Success rate: {stats['success_rate_overall']:.1f}%")

print("\n🧪 Step 3: Get recommendations for your reaction")
print("   Example: Bromobenzene + Ethylamine (C-N Coupling)")
print("   Code: result = recommender.recommend('c1ccc(Br)cc1', 'CCN', top_k=3)")

result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles="CCN",
    top_k=3
)

print(f"\n✅ Results:")
print(f"   Reactant types detected: {result.reactant_a_type} + {result.reactant_b_type}")
print(f"   Predicted reaction: {result.predicted_reaction_type}")
print(f"   Matching experiments: {result.total_matching_experiments}")
print(f"   Recommendations found: {len(result.recommendations)}")

if result.recommendations:
    rec = result.recommendations[0]
    print(f"\n🏆 Top Recommendation:")
    print(f"   Confidence: {rec.confidence_score:.1f}/100")
    print(f"   Catalyst: {rec.catalyst}")
    print(f"   Ligand: {rec.ligand}")
    print(f"   Base: {rec.base}")
    print(f"   Solvent: {rec.solvent}")
    print(f"   Success rate: {rec.success_rate:.1f}% ({rec.num_experiments} experiments)")
    print(f"   Average yield: {rec.avg_yield:.1f}%")

print("\n" + "="*80)
print("💡 USAGE TIPS")
print("="*80)
print("""
Python API:
  from chemtools.HTE import HTERecommender
  recommender = HTERecommender()
  result = recommender.recommend(reactant_a, reactant_b, top_k=5)

Command Line:
  # Simple query
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN"
  
  # Compact output
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN" --compact
  
  # JSON format
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN" --json
  
  # Filter by reaction type
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN" --reaction C_N_Coupling
  
  # Batch processing
  python -m chemtools.HTE.cli --batch examples/hte_queries.txt

Batch File Format (examples/hte_queries.txt):
  c1ccc(Br)cc1 CCN
  c1ccc(Cl)cc1 c1ccc(N)cc1
  c1ccc(Br)cc1 c1ccc(B(O)O)cc1
""")

print("\n📚 See docs/HTE_RECOMMENDER.md for complete documentation")
print("="*80)
