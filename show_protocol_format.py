"""
Show the complete protocol recommendation output format
"""

import json
from chemtools.protocol.recommend import ProtocolRecommender

# Initialize recommender
recommender = ProtocolRecommender()

# Test reaction
reaction = "O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1"

# Get recommendation
results = recommender.recommend(
    reaction_smiles=reaction,
    k=1,
    use_smarts_filter=True,
    use_standard_format=True
)

# Pretty print
print("=" * 80)
print("PROTOCOL RECOMMENDATION OUTPUT (Standard Format)")
print("=" * 80)
print(json.dumps(results, indent=2))

# Show just one recommendation entry in detail
if results.get('recommended_conditions'):
    print("\n" + "=" * 80)
    print("SINGLE RECOMMENDATION ENTRY (Detailed)")
    print("=" * 80)
    rec = results['recommended_conditions'][0]
    print(json.dumps(rec, indent=2))
