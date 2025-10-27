"""Fix incorrect entity_id values in aliases.json."""

import json

# Mapping from incorrect to correct entity IDs
FIXES = {
    'alkene': 'Alkene',
    'allylic_halide': 'Allylic-X',
    'benzylic_halide': 'Benzylic-X',
    'imines': 'Imine',
    'acyl_source_electrophile': 'Acyl-electrophile',
    'rso2cl': 'RSO2Cl',
}

# Load aliases
with open('chemtools/taxonomy/data/aliases.json', 'r', encoding='utf-8') as f:
    aliases = json.load(f)

print(f"Loaded {len(aliases)} aliases")
print()

# Fix entity_ids
fixes_applied = 0
for alias_obj in aliases:
    entity_id = alias_obj.get('entity_id')
    if entity_id in FIXES:
        old_id = entity_id
        new_id = FIXES[entity_id]
        alias_obj['entity_id'] = new_id
        fixes_applied += 1
        print(f"Fixed: {alias_obj.get('alias'):30s} {old_id:30s} → {new_id}")

print()
print(f"Total fixes applied: {fixes_applied}")

# Also check for the encoding issue
encoding_issues = []
for alias_obj in aliases:
    alias_name = alias_obj.get('alias', '')
    entity_id = alias_obj.get('entity_id', '')
    if '酶' in alias_name or '酶' in entity_id:  # The mojibake character
        encoding_issues.append(alias_obj)

if encoding_issues:
    print(f"\nFound {len(encoding_issues)} aliases with encoding issues:")
    for a in encoding_issues:
        print(f"  Alias: {a.get('alias'):40s} → Entity: {a.get('entity_id')}")
    
    # Fix the specific ones we know about
    for alias_obj in encoding_issues:
        if 'phosphoric_acid' in alias_obj.get('alias', '') or 'phosphoric_acid' in alias_obj.get('entity_id', ''):
            old_alias = alias_obj.get('alias')
            old_entity = alias_obj.get('entity_id')
            # The correct term is "brønsted" (with ø), but use ASCII-safe "bronsted"
            alias_obj['alias'] = 'phosphoric_acid_bronsted'
            alias_obj['entity_id'] = 'phosphoric_acid_bronsted'
            print(f"\nFixed encoding:")
            print(f"  Alias: {old_alias} → {alias_obj['alias']}")
            print(f"  Entity: {old_entity} → {alias_obj['entity_id']}")
            fixes_applied += 1

# Save the fixed file
with open('chemtools/taxonomy/data/aliases.json', 'w', encoding='utf-8') as f:
    json.dump(aliases, f, indent=4, ensure_ascii=False)

print(f"\n✓ Saved {len(aliases)} aliases with {fixes_applied} fixes")
