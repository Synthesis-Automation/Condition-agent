import json

with open('chemtools/taxonomy/data/reactant_types.json', encoding='utf-8') as f:
    data = json.load(f)

print(f'Total reactant type categories: {len(data)}')
print('\nExisting categories:')
for i, item in enumerate(data):
    name = item.get('name', 'N/A')
    id_val = item.get('id', 'N/A')
    num_members = len(item.get('members', []))
    print(f'{i+1:2d}. {id_val:35s} | {name:25s} | {num_members} members')
