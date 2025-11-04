"""Validate the expanded calculable_features.json"""
import json
from pathlib import Path

json_file = Path(__file__).parent.parent / "chemtools" / "featurizers" / "calculable_features.json"

with open(json_file, 'r', encoding='utf-8') as f:
    data = json.load(f)

print("✓ Valid JSON!")
print(f"\nVersion: {data['version']}")
print(f"Description: {data.get('description', 'N/A')[:80]}...")
print(f"\nFeatures: {len(data['features'])}")
print(f"Derived shortcuts: {len(data['derived_shortcuts'])}")
print(f"Categories: {len(data['schema_notes']['categories'])}")

print("\n=== Feature Categories ===")
for cat in data['schema_notes']['categories']:
    print(f"  - {cat}")

print("\n=== Sample Features (first 10) ===")
for i, feat in enumerate(data['features'][:10], 1):
    category = feat.get('category', 'N/A')
    print(f"{i}. {feat['token']} ({feat['type']}) - {category}")
    print(f"   Why: {feat.get('why', 'N/A')[:60]}...")

print("\n=== Sample Derived Shortcuts (first 5) ===")
for i, derived in enumerate(data['derived_shortcuts'][:5], 1):
    print(f"{i}. {derived['token']}")
    print(f"   Derive: {derived.get('derive', 'N/A')[:80]}...")
    print(f"   Why: {derived.get('why', 'N/A')[:60]}...")

print("\n=== Category Distribution ===")
from collections import Counter
category_counts = Counter(feat.get('category', 'uncategorized') for feat in data['features'])
for cat, count in sorted(category_counts.items(), key=lambda x: -x[1])[:15]:
    print(f"  {cat}: {count} features")

print(f"\n✓ Expansion successful!")
print(f"  Total coverage: {len(data['features']) + len(data['derived_shortcuts'])} tokens")
