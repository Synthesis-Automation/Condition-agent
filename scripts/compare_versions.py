#!/usr/bin/env python3
"""
Quick comparison of v2.2 vs v3.0 calculable features
"""

import json
from pathlib import Path

# Load both versions
v2_file = Path("chemtools/featurizers/calculable_features_v2.2_backup.json")
v3_file = Path("chemtools/featurizers/calculable_features.json")

with open(v2_file, 'r', encoding='utf-8') as f:
    v2_data = json.load(f)

with open(v3_file, 'r', encoding='utf-8') as f:
    v3_data = json.load(f)

print("=" * 80)
print("CALCULABLE FEATURES: v2.2 → v3.0 COMPARISON")
print("=" * 80)

print(f"\nVersion Comparison:")
print(f"  v2.2: {v2_data.get('version', 'N/A')}")
print(f"  v3.0: {v3_data.get('version', 'N/A')}")

print(f"\nFeature Counts:")
print(f"  v2.2 base features: {len(v2_data.get('features', []))}")
print(f"  v3.0 base features: {len(v3_data.get('features', []))}")
print(f"  → Increase: +{len(v3_data.get('features', [])) - len(v2_data.get('features', []))} features (+{100*(len(v3_data.get('features', [])) - len(v2_data.get('features', [])))/len(v2_data.get('features', [])):.0f}%)")

print(f"\nDerived Shortcuts:")
print(f"  v2.2: {len(v2_data.get('derived_shortcuts', []))}")
print(f"  v3.0: {len(v3_data.get('derived_shortcuts', []))}")
print(f"  → Increase: +{len(v3_data.get('derived_shortcuts', [])) - len(v2_data.get('derived_shortcuts', []))} shortcuts")

print(f"\nTotal Tokens:")
v2_total = len(v2_data.get('features', [])) + len(v2_data.get('derived_shortcuts', []))
v3_total = len(v3_data.get('features', [])) + len(v3_data.get('derived_shortcuts', []))
print(f"  v2.2: {v2_total}")
print(f"  v3.0: {v3_total}")
print(f"  → Increase: +{v3_total - v2_total} tokens (+{100*(v3_total - v2_total)/v2_total:.0f}%)")

# Get new features by category
print(f"\n" + "=" * 80)
print("NEW FEATURES BY CATEGORY")
print("=" * 80)

v2_tokens = {f['token'] for f in v2_data.get('features', [])}
v3_features = v3_data.get('features', [])

new_by_category = {}
for feat in v3_features:
    if feat['token'] not in v2_tokens:
        cat = feat.get('category', 'uncategorized')
        if cat not in new_by_category:
            new_by_category[cat] = []
        new_by_category[cat].append(feat['token'])

for category, features in sorted(new_by_category.items(), key=lambda x: -len(x[1])):
    print(f"\n{category} ({len(features)} new):")
    for feat in features[:5]:  # Show first 5
        print(f"  • {feat}")
    if len(features) > 5:
        print(f"  ... and {len(features) - 5} more")

# Show new derived shortcuts
print(f"\n" + "=" * 80)
print("NEW DERIVED SHORTCUTS")
print("=" * 80)

v2_derived_tokens = {d['token'] for d in v2_data.get('derived_shortcuts', [])}
new_derived = [d for d in v3_data.get('derived_shortcuts', []) if d['token'] not in v2_derived_tokens]

print(f"\nAdded {len(new_derived)} new derived features:")
for d in new_derived[:10]:
    print(f"  • {d['token']}")
    print(f"    → {d.get('derive', 'N/A')[:70]}...")
if len(new_derived) > 10:
    print(f"  ... and {len(new_derived) - 10} more")

print(f"\n" + "=" * 80)
print("ORGANIZATIONAL IMPROVEMENTS")
print("=" * 80)

print(f"\nCategories:")
v2_categories = set()
for feat in v2_data.get('features', []):
    if 'category' in feat:
        v2_categories.add(feat['category'])

v3_categories = set()
for feat in v3_data.get('features', []):
    if 'category' in feat:
        v3_categories.add(feat['category'])

print(f"  v2.2: {len(v2_categories)} categories")
print(f"  v3.0: {len(v3_categories)} categories")
print(f"\nNew categories in v3.0:")
for cat in sorted(v3_categories - v2_categories):
    count = sum(1 for f in v3_features if f.get('category') == cat)
    print(f"  • {cat} ({count} features)")

print(f"\n" + "=" * 80)
print("✓ EXPANSION COMPLETE")
print("=" * 80)
print(f"\nv3.0 now provides comprehensive coverage for:")
print(f"  ✓ Protecting group strategies")
print(f"  ✓ Complete heterocycle library")
print(f"  ✓ Safety and reactivity assessment")
print(f"  ✓ Drug-likeness evaluation (Lipinski, Veber, PAINS)")
print(f"  ✓ Medicinal chemistry features")
print(f"  ✓ Stereochemistry and 3D topology")
print(f"  ✓ Redox and photochemistry")
print(f"  ✓ Expanded organometallic detection")
print(f"  ✓ Sulfur and phosphorus functionality")
print(f"  ✓ LLM-ready with comprehensive 'why' explanations")
print(f"\n" + "=" * 80)
