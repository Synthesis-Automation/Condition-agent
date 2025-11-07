#!/usr/bin/env python3
"""Check calculable_features.json for logical issues and redundancies."""

import json
from collections import defaultdict

# Load the JSON file
with open('chemtools/featurizers/calculable_features.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

features = data.get('features', [])
derived = data.get('derived_shortcuts', [])

print(f"Total features: {len(features)}")
print(f"Derived shortcuts: {len(derived)}")
print("\n" + "="*80)

# 1. Check for duplicate tokens
tokens = [f.get('token') for f in features]
duplicates = [t for t in set(tokens) if tokens.count(t) > 1]
if duplicates:
    print(f"\n❌ DUPLICATE TOKENS FOUND ({len(duplicates)}):")
    for t in duplicates:
        print(f"  - {t} (appears {tokens.count(t)} times)")
else:
    print("\n✅ No duplicate tokens")

# 2. Check for tokens without detection methods
no_detection = [f for f in features if not f.get('detect')]
if no_detection:
    print(f"\n❌ FEATURES WITHOUT DETECTION METHOD ({len(no_detection)}):")
    for f in no_detection[:10]:
        print(f"  - {f.get('token')}")
else:
    print("\n✅ All features have detection methods")

# 3. Analyze feature categories
categories = defaultdict(list)
for feat in features:
    tok = feat.get('token', '')
    if 'halide' in tok or 'chloride' in tok or 'bromide' in tok or 'iodide' in tok or 'fluoride' in tok:
        categories['halides'].append(tok)
    elif 'amine' in tok or 'amino' in tok:
        categories['amines'].append(tok)
    elif 'alcohol' in tok or 'hydroxyl' in tok:
        categories['alcohols'].append(tok)
    elif 'carbonyl' in tok or 'ketone' in tok or 'aldehyde' in tok:
        categories['carbonyls'].append(tok)
    elif 'ester' in tok:
        categories['esters'].append(tok)
    elif 'acid' in tok:
        categories['acids'].append(tok)
    elif 'protecting' in tok or 'protected' in tok or 'boc' in tok or 'cbz' in tok or 'fmoc' in tok:
        categories['protecting_groups'].append(tok)
    elif 'boronic' in tok or 'boron' in tok:
        categories['boron'].append(tok)

print("\n" + "="*80)
print("\nFEATURE CATEGORIES:")
for cat, feats in sorted(categories.items()):
    print(f"\n{cat.upper()} ({len(feats)}):")
    for f in sorted(feats):
        print(f"  - {f}")

# 4. Check for potential redundancies
print("\n" + "="*80)
print("\nPOTENTIAL REDUNDANCIES:")

# Check if aryl_halide_present could be redundant with sp2 halides
aryl_halide = next((f for f in features if f.get('token') == 'aryl_halide_present'), None)
sp2_halides = [f for f in features if f.get('token', '').startswith('sp2_') and 'halide' in f.get('token', '')]

if aryl_halide and sp2_halides:
    print(f"\n⚠️  aryl_halide_present vs sp2_*_present:")
    print(f"  - aryl_halide_present: {aryl_halide.get('why', '')}")
    print(f"  - sp2 halides: {[f.get('token') for f in sp2_halides]}")
    print("  → These could overlap - aryl halides are a subset of sp2 halides")

# Check amine hierarchy
amine_features = [f for f in features if 'amine' in f.get('token', '')]
if amine_features:
    print(f"\n⚠️  Amine hierarchy ({len(amine_features)} features):")
    for f in sorted(amine_features, key=lambda x: x.get('token')):
        print(f"  - {f.get('token')}: {f.get('why', '')}")

# 5. Check SMARTS patterns for issues
print("\n" + "="*80)
print("\nSMART PATTERN ANALYSIS:")

carboxylic_acid = next((f for f in features if f.get('token') == 'carboxylic_acid_present'), None)
if carboxylic_acid:
    smarts = carboxylic_acid.get('detect', {}).get('smarts_any', [])
    print(f"\ncarboxylic_acid_present SMARTS: {smarts}")
    if 'C(=O)O[H]' in smarts:
        print("  ⚠️  WARNING: Pattern C(=O)O[H] requires explicit H - may miss standard SMILES")
    if 'C(=O)[OH]' in smarts:
        print("  ✅ Pattern C(=O)[OH] is correct - matches with or without explicit H")

# 6. Check for missing common features
print("\n" + "="*80)
print("\nCOMMON FEATURE CHECK:")

common_features = [
    'primary_amine_present',
    'secondary_amine_present', 
    'tertiary_amine_present',
    'primary_alcohol_present',
    'secondary_alcohol_present',
    'tertiary_alcohol_present',
    'carboxylic_acid_present',
    'ester_present',
    'amide_present',
    'alkene_present',
    'alkyne_present',
    'aromatic_present'
]

for expected in common_features:
    found = any(f.get('token') == expected for f in features)
    status = "✅" if found else "❌"
    print(f"{status} {expected}")

print("\n" + "="*80)
print("\nANALYSIS COMPLETE")
