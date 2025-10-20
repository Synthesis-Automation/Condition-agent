import json

with open('protocol_recommender_test_results.json', 'r') as f:
    data = json.load(f)

# Find ALL Sonogashira reactions
sonogashira = [r for r in data['results'] if r['category'] == 'Sonogashira']

print("=" * 80)
print(f"ALL SONOGASHIRA REACTIONS ({len(sonogashira)} total)")
print("=" * 80)
print()

for i, r in enumerate(sonogashira, 1):
    status_emoji = "✓" if r['status'] == 'matched' else "✗"
    print(f"{status_emoji} {i}. [{r['status']:9s}] {r['description'][:60]}")
    print(f"   SMILES: {r['smiles']}")
    if r['status'] == 'matched':
        print(f"   Top match: {r['top_family']} (sim={r['top_similarity']:.3f})")
    print()
