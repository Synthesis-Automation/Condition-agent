"""Find and fix problematic aliases in reactant_types.json."""

import json

# Load the file
with open('chemtools/taxonomy/data/reactant_types.json', encoding='utf-8') as f:
    data = json.load(f)

# The file is a list of reactant type objects
reactant_types = data if isinstance(data, list) else data.get('reactant_types', [])

# Get valid IDs
valid_ids = {rt['id'] for rt in reactant_types if 'id' in rt}
print(f"Valid reactant type IDs: {len(valid_ids)}")
print()

# Mapping from incorrect aliases to correct IDs
ALIAS_FIXES = {
    'alkene': 'Alkene',
    'allylic_halide': 'Allylic-X',
    'benzylic_halide': 'Benzylic-X',
    'imines': 'Imine',
    'acyl_source_electrophile': 'Acyl-electrophile',
    'rso2cl': 'RSO2Cl',
}

# Find all aliases in the file
all_aliases = []
for rt in reactant_types:
    for alias_obj in rt.get('aliases', []):
        alias_obj['parent_id'] = rt['id']
        all_aliases.append(alias_obj)
    for member in rt.get('members', []):
        for alias_obj in member.get('aliases', []):
            alias_obj['parent_id'] = member['id']
            all_aliases.append(alias_obj)

print(f"Total aliases found: {len(all_aliases)}")
print()

# Find problematic aliases
problematic = []
for alias_obj in all_aliases:
    target = alias_obj.get('target_id')
    if target and target in ALIAS_FIXES:
        problematic.append(alias_obj)

print(f"Found {len(problematic)} aliases with incorrect targets")
for a in problematic[:30]:
    old_target = a.get('target_id', 'N/A')
    new_target = ALIAS_FIXES.get(old_target, old_target)
    alias_name = a.get('alias', a.get('name', 'N/A'))
    print(f"  {alias_name:30s} -> {old_target:30s} (should be: {new_target})")
