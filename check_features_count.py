import json

with open('chemtools/featurizers/calculable_features.json') as f:
    data = json.load(f)

print(f"Features: {len(data['features'])}")
print(f"Derived: {len(data.get('derived', []))}")
print(f"Version: {data.get('version')}")

# Count reactant features
reactant_count = sum(1 for f in data['features'] if f['token'].endswith('_reactant'))
present_count = sum(1 for f in data['features'] if f['token'].endswith('_present'))

print(f"\nReactant features (_reactant suffix): {reactant_count}")
print(f"Present features (_present suffix): {present_count}")
print(f"Other features: {len(data['features']) - reactant_count - present_count}")
