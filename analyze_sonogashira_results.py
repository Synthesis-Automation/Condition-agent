import json

with open('protocol_recommender_test_results.json', 'r') as f:
    data = json.load(f)

# Find matched Sonogashira reactions
matches = [r for r in data['results'] if r['category'] == 'Sonogashira' and r['status'] == 'matched']
no_matches = [r for r in data['results'] if r['category'] == 'Sonogashira' and r['status'] == 'no_match']

print("=" * 80)
print(f"SONOGASHIRA REACTIONS ANALYSIS")
print("=" * 80)
print()

print(f"Matched Sonogashira reactions ({len(matches)} of 8):")
print("-" * 80)
for i, r in enumerate(matches, 1):
    print(f"{i}. {r['description']}")
    print(f"   SMILES: {r['smiles']}")
    print(f"   Top match: {r['top_family']} (similarity={r['top_similarity']:.3f})")
    print()

print("=" * 80)
print(f"Un-matched Sonogashira reactions ({len(no_matches)} of 8):")
print("-" * 80)
for i, r in enumerate(no_matches, 1):
    print(f"{i}. {r['description']}")
    print(f"   SMILES: {r['smiles']}")
    print()
