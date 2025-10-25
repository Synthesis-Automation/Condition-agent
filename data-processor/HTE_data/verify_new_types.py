import json

# Load the updated reactant_types.json
with open('data-processor/other_data/reactant_types.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

print(f'Total categories: {len(data)}')
total_members = sum(len(v['members']) for v in data.values())
print(f'Total member types: {total_members}')

print('\n=== New Categories Added ===')
new_cats = ['Heterocyclic-halide', 'Enamine', 'Imines', 'Allylic-halide', 'Benzylic-halide', 'Azide', 'Nitrile']
for cat in new_cats:
    if cat in data:
        print(f'\n{cat}:')
        print(f'  Members: {len(data[cat]["members"])}')
        print(f'  Group: {data[cat]["group"]}')
        print(f'  Description: {data[cat]["description"]}')
        for member in data[cat]['members']:
            print(f'    - {member["id"]}: {member["name"]}')
    else:
        print(f'\n{cat}: NOT FOUND')

print('\n=== All Categories by Group ===')
groups = {}
for cat, info in data.items():
    group = info['group']
    if group not in groups:
        groups[group] = []
    groups[group].append(cat)

for group, cats in sorted(groups.items()):
    print(f'\n{group}:')
    for cat in sorted(cats):
        member_count = len(data[cat]['members'])
        print(f'  - {cat} ({member_count} members)')
