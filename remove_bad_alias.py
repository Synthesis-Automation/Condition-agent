"""Remove the problematic phosphoric_acid alias."""

import json

# Load aliases
with open('chemtools/taxonomy/data/aliases.json', 'r', encoding='utf-8') as f:
    aliases = json.load(f)

print(f"Loaded {len(aliases)} aliases")

# Find and remove problematic alias
before = len(aliases)
aliases = [a for a in aliases if 'phosphoric_acid' not in a.get('alias', '').lower()]
after = len(aliases)

removed = before - after
print(f"Removed {removed} aliases containing 'phosphoric_acid'")

# Save
with open('chemtools/taxonomy/data/aliases.json', 'w', encoding='utf-8') as f:
    json.dump(aliases, f, indent=4, ensure_ascii=False)

print(f"✓ Saved {len(aliases)} aliases")
