from chemtools.protocol.recommend import ProtocolRecommender
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

# Test the iodide case
reaction_smiles = "Ic1ccncc1.C#CC>>C#Cc1ccncc1"

print("=" * 80)
print("TESTING IODIDE SONOGASHIRA WITH FULL DEBUG")
print("=" * 80)
print()
print(f"Reaction: {reaction_smiles}")
print()

recommender = ProtocolRecommender()

# Test with SMARTS filtering
print("Testing with SMARTS filter enabled...")
result = recommender.recommend(reaction_smiles, k=5, use_smarts_filter=True)

print()
# Standard format uses 'recommended_conditions', not 'protocols'
recommendations = result.get('recommended_conditions', [])
print(f"Number of matches: {len(recommendations)}")
print()

if recommendations:
    for i, rec in enumerate(recommendations[:3], 1):
        metadata = rec.get('protocol_metadata', {})
        print(f"{i}. {metadata.get('reaction_family', 'Unknown')}")
        print(f"   Confidence: {rec.get('confidence', 0):.3f}")
        print(f"   File: {metadata.get('filename', 'N/A')}")
        print()
else:
    print("No matches found!")
    if 'extras' in result:
        print(f"Extras: {result['extras']}")

print()
print("=" * 80)
