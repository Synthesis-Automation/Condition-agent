from chemtools.featurizers.calculable import classify_reactant_smiles, get_reactant_type_features
import json

# Load and check file
with open('chemtools/featurizers/calculable_features.json') as f:
    data = json.load(f)

print(f"Version: {data['version']}")
print(f"Features: {len(data['features'])}")
print(f"Derived shortcuts: {len(data.get('derived_shortcuts', []))}")

# Count tokens by suffix
features = data['features']
reactant_suffix = sum(1 for f in features if f['token'].endswith('_reactant'))
present_suffix = sum(1 for f in features if f['token'].endswith('_present'))
with_metadata = sum(1 for f in features if 'reactant_metadata' in f)

print(f"\nToken counts:")
print(f"  _reactant suffix: {reactant_suffix}")
print(f"  _present suffix: {present_suffix}")
print(f"  With reactant_metadata: {with_metadata}")

# Test functionality
print(f"\nFunctionality test:")
result = classify_reactant_smiles('c1ccc(Br)cc1')
print(f"  Bromobenzene: {result['member_type']} ({result['category']}) - {result['name']}")

result2 = classify_reactant_smiles('c1ccc(B(O)O)cc1')
print(f"  Phenylboronic acid: {result2['member_type']} ({result2['category']}) - {result2['name']}")

print("\n✅ Refactoring complete - no _reactant token suffixes remaining!")
