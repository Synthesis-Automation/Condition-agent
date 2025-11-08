"""Quick check of feature counts in calculable_features.json"""
import json

with open('chemtools/featurizers/calculable_features.json', 'r', encoding='utf-8') as f:
    spec = json.load(f)

print(f"Version: {spec['version']}")
print(f"Total base features: {len(spec['features'])}")
derived = spec.get('derived_shortcuts') or []
print(f"Derived features: {len(derived)}")
print(f"Grand total: {len(spec['features']) + len(derived)}")

# Count by type
bool_count = sum(1 for f in spec['features'] if f.get('type') == 'bool')
int_count = sum(1 for f in spec['features'] if f.get('type') == 'int')

print(f"\nBreakdown:")
print(f"  Boolean features: {bool_count}")
print(f"  Integer features: {int_count}")
print(f"  Derived features: {len(derived)}")
