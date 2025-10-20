from chemtools.protocol.recommend import ProtocolRecommender
import logging

# Enable info logging
logging.basicConfig(level=logging.INFO)

# Test bromide Sonogashira
reaction_smiles = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"

print("=" * 80)
print("TESTING BROMIDE SONOGASHIRA (Diphenylacetylene)")
print("=" * 80)
print()
print(f"Reaction: {reaction_smiles}")
print()

recommender = ProtocolRecommender()

# Test with SMARTS filtering
print("Testing with SMARTS filter enabled...")
result = recommender.recommend(reaction_smiles, k=5, use_smarts_filter=True)

print()
recommendations = result.get('recommended_conditions', [])
print(f"Number of matches: {len(recommendations)}")
print()

if recommendations:
    for i, rec in enumerate(recommendations[:5], 1):
        metadata = rec.get('protocol_metadata', {})
        print(f"{i}. {metadata.get('reaction_family', 'Unknown')}")
        print(f"   Confidence: {rec.get('confidence', 0):.3f}")
        print(f"   File: {metadata.get('filename', 'N/A')}")
        print(f"   SMARTS: {metadata.get('reaction_smarts', [])}")
        print()
else:
    print("No matches found!")
    if 'extras' in result:
        print(f"Extras: {result['extras']}")

print()
print("=" * 80)
