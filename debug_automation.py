"""
Direct test to debug automation format
"""

from chemtools.recommend import UnifiedRecommender

recommender = UnifiedRecommender()

# Sonogashira reaction
rxn = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"

print("Testing with format_for_automation=True, source_types=['rule']")
print()

results = recommender.recommend(
    reaction_smiles=rxn,
    top_k=3,
    source_types=['rule'],
    format_for_automation=True,
    scale_mmol=2.0
)

print(f"Found {len(results)} results")
print()

for r in results:
    print(f"[{r.rank}] {r.name}")
    print(f"  Type: {r.source_type}")
    print(f"  Similarity: {r.similarity:.3f}")
    print(f"  Has full_data attribute: {hasattr(r, 'full_data')}")
    
    if hasattr(r, 'full_data'):
        print(f"  full_data is not None: {r.full_data is not None}")
        if r.full_data:
            print(f"  full_data keys: {list(r.full_data.keys())}")
            if 'reaction_setup' in r.full_data:
                setup = r.full_data['reaction_setup'][0]
                chems = setup.get('chemicals', [])
                print(f"  ✅ Has reaction_setup with {len(chems)} chemicals")
                for i, chem in enumerate(chems[:3], 1):
                    print(f"      {i}. {chem['name']} ({chem['role']})")
            else:
                print(f"  ❌ No reaction_setup in full_data")
        else:
            print(f"  ❌ full_data is None")
    else:
        print(f"  ❌ No full_data attribute")
    print()
