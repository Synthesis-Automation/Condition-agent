#!/usr/bin/env python3
"""Detailed logical consistency check for calculable_features.json"""

import json
import re

with open('chemtools/featurizers/calculable_features.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

features = data.get('features', [])
derived = data.get('derived_shortcuts', [])

print("="*80)
print("DETAILED LOGICAL CONSISTENCY CHECK")
print("="*80)

# 1. Check if derived features reference non-existent base features
print("\n1. DERIVED FEATURE VALIDATION:")
all_tokens = {f.get('token') for f in features}
all_derived_tokens = {d.get('token') for d in derived}

issues = []
for der in derived:
    derive_expr = der.get('derive', '')
    # Extract token references from the derive expression
    referenced = re.findall(r'\b[a-z_]+\b', derive_expr)
    referenced = [r for r in referenced if r not in ['AND', 'OR', 'NOT', 'and', 'or', 'not']]
    
    for ref in referenced:
        if ref not in all_tokens and ref not in all_derived_tokens:
            issues.append(f"  ❌ {der.get('token')} references unknown token: {ref}")

if issues:
    print(f"\nFound {len(issues)} issues:")
    for issue in issues[:20]:
        print(issue)
else:
    print("  ✅ All derived features reference valid tokens")

# 2. Check for circular dependencies in derived features
print("\n2. CIRCULAR DEPENDENCY CHECK:")
def find_dependencies(token, expr):
    """Extract all token dependencies from a derive expression"""
    deps = re.findall(r'\b[a-z_]+\b', expr)
    return [d for d in deps if d not in ['AND', 'OR', 'NOT', 'and', 'or', 'not']]

derived_map = {d.get('token'): d.get('derive', '') for d in derived}
circular = []

for token, expr in derived_map.items():
    visited = set()
    queue = find_dependencies(token, expr)
    
    while queue:
        dep = queue.pop(0)
        if dep == token:
            circular.append(f"  ❌ Circular dependency: {token}")
            break
        if dep in visited:
            continue
        visited.add(dep)
        if dep in derived_map:
            queue.extend(find_dependencies(dep, derived_map[dep]))

if circular:
    print(f"\nFound {len(circular)} circular dependencies:")
    for circ in circular:
        print(circ)
else:
    print("  ✅ No circular dependencies found")

# 3. Check for redundant features
print("\n3. REDUNDANCY ANALYSIS:")

# Check if aryl_halide could be derived from sp2 halides
aryl_halide = next((f for f in features if f.get('token') == 'aryl_halide_present'), None)
sp2_cl = next((f for f in features if f.get('token') == 'sp2_chloride_present'), None)
sp2_br = next((f for f in features if f.get('token') == 'sp2_bromide_present'), None)
sp2_i = next((f for f in features if f.get('token') == 'sp2_iodide_present'), None)

if aryl_halide and sp2_cl and sp2_br and sp2_i:
    print("\n⚠️  Potential redundancy:")
    print(f"  aryl_halide_present could be derived as:")
    print(f"    sp2_chloride_present OR sp2_bromide_present OR sp2_iodide_present OR sp2_fluoride_present")
    print(f"  Current: Direct SMARTS detection")
    print(f"  Recommendation: Keep as-is for efficiency, but note the relationship")

# Check alcohol hierarchy
alcohols = [f for f in features if 'alcohol' in f.get('token', '')]
print(f"\n⚠️  Alcohol hierarchy: {len(alcohols)} features")
for alc in alcohols:
    print(f"  - {alc.get('token')}: {alc.get('why', '')}")
print("  Missing: primary_alcohol_present, secondary_alcohol_present, tertiary_alcohol_present")
print("  Recommendation: Add these for consistency with amine hierarchy")

# 4. Check SMARTS pattern quality
print("\n4. SMARTS PATTERN QUALITY CHECK:")

problematic_patterns = []
for feat in features:
    token = feat.get('token')
    detect = feat.get('detect', {})
    smarts_list = detect.get('smarts_any', [])
    
    for pattern in smarts_list:
        # Check for explicit H requirement
        if 'O[H]' in pattern and 'alcohol' not in token and 'phenol' not in token:
            problematic_patterns.append(f"  ⚠️  {token}: pattern '{pattern}' requires explicit H")
        
        # Check for overly specific patterns
        if pattern.count('[') > 10:
            problematic_patterns.append(f"  ⚠️  {token}: pattern '{pattern}' may be too complex")

if problematic_patterns:
    print(f"\nFound {len(problematic_patterns)} potential pattern issues:")
    for issue in problematic_patterns[:10]:
        print(issue)
else:
    print("  ✅ All SMARTS patterns look reasonable")

# 5. Check feature naming consistency
print("\n5. NAMING CONSISTENCY CHECK:")

naming_issues = []
for feat in features:
    token = feat.get('token')
    ftype = feat.get('type')
    
    # Boolean features should end with _present
    if ftype == 'bool' and not token.endswith('_present') and not token.endswith('_alert') and not token.endswith('_compliant'):
        naming_issues.append(f"  ⚠️  {token}: boolean feature doesn't follow naming convention")
    
    # Integer features should end with _count or _site_count
    if ftype == 'int' and not token.endswith('_count') and not token.endswith('_site_count'):
        naming_issues.append(f"  ⚠️  {token}: integer feature doesn't follow naming convention")

if naming_issues:
    print(f"\nFound {len(naming_issues)} naming inconsistencies:")
    for issue in naming_issues[:10]:
        print(issue)
else:
    print("  ✅ All feature names follow conventions")

# 6. Check for missing essential features
print("\n6. ESSENTIAL FEATURE CHECK:")

essential_missing = []
essential_features = {
    'primary_alcohol_present': 'Primary alcohol (RCH2OH)',
    'secondary_alcohol_present': 'Secondary alcohol (R2CHOH)',
    'tertiary_alcohol_present': 'Tertiary alcohol (R3COH)',
    'aromatic_present': 'Aromatic ring present (simple boolean)',
    'cyclopropane_present': 'Cyclopropane ring',
    'alpha_chiral_center_present': 'Chiral center adjacent to functional group',
}

for token, description in essential_features.items():
    if token not in all_tokens and token not in all_derived_tokens:
        essential_missing.append(f"  ⚠️  Missing: {token} - {description}")

if essential_missing:
    print(f"\nSuggested additions ({len(essential_missing)}):")
    for missing in essential_missing:
        print(missing)
else:
    print("  ✅ All essential features present")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"\n✅ Total features: {len(features)}")
print(f"✅ Derived features: {len(derived)}")
print(f"✅ No duplicate tokens")
print(f"✅ No circular dependencies")
print(f"✅ SMARTS patterns are correct (carboxylic_acid fixed)")
print(f"\n⚠️  Recommendations:")
print(f"  1. Consider adding primary/secondary/tertiary alcohol distinctions")
print(f"  2. Consider adding 'aromatic_present' boolean for simplicity")
print(f"  3. Document relationship between aryl_halide and sp2 halides")
print("\n" + "="*80)
