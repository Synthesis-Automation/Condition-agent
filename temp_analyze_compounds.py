import json

# Read the JSON file
with open(r'c:\Git-softwares\Condition-agent\chemtools\taxonomy\data\organic_compounds.v1.3.json', encoding='utf-8') as f:
    data = json.load(f)

# Find all mismatches where id != name
mismatches = []
for compound in data['compounds']:
    if compound['id'] != compound['name']:
        mismatches.append({
            'id': compound['id'],
            'name': compound['name'],
            'A': compound.get('A', 'N/A'),
            'B': compound.get('B', 'N/A')
        })

# Print results
print(f"Found {len(mismatches)} compounds where name != id\n")
print("="*70)

for m in mismatches:
    print(f"ID:   {m['id']}")
    print(f"Name: {m['name']}")
    print(f"A:    {m['A']}")
    print(f"B:    {m['B']}")
    print("-"*70)

print(f"\nTotal compounds analyzed: {len(data['compounds'])}")
print(f"Compounds with matching id/name: {len(data['compounds']) - len(mismatches)}")
print(f"Compounds with mismatched id/name: {len(mismatches)}")
