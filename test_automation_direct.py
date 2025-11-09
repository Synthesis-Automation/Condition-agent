"""
Simple direct test of automation format.
"""

from chemtools.recommend import UnifiedRecommender

print("=" * 80)
print("DIRECT TEST: Automation Format")
print("=" * 80)
print()

# Initialize
recommender = UnifiedRecommender()

# Test reaction
rxn = "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"

print(f"Reaction: {rxn}")
print()

# Get with automation format
print("Requesting with format_for_automation=True, scale_mmol=2.0...")
results = recommender.recommend(
    reaction_smiles=rxn,
    top_k=2,
    format_for_automation=True,
    scale_mmol=2.0
)

print(f"Found {len(results)} results\n")

for result in results:
    print(f"[{result.rank}] {result.name}")
    print(f"    Similarity: {result.similarity:.3f}")
    print(f"    Type: {result.source_type}")
    
    if hasattr(result, 'full_data') and result.full_data:
        print(f"    ✅ Has automation format data")
        
        if 'reaction_setup' in result.full_data:
            setup = result.full_data['reaction_setup'][0]
            chemicals = setup.get('chemicals', [])
            print(f"    Addition sequence ({len(chemicals)} steps):")
            for i, chem in enumerate(chemicals[:5], 1):  # Show first 5
                name = chem['name']
                role = chem['role']
                print(f"      {i}. {name} ({role})")
            if len(chemicals) > 5:
                print(f"      ... +{len(chemicals)-5} more")
    else:
        print(f"    ❌ No automation format data")
    
    print()

print("=" * 80)
print("✅ TEST COMPLETE")
print("=" * 80)
