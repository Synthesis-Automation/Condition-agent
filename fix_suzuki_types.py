"""Quick fix for Suzuki_db.json mixed types."""

import json
from pathlib import Path

rule_file = Path('data/rule_db_v2/Suzuki_db.json')

print("=" * 80)
print("FIXING SUZUKI_DB.JSON MIXED TYPES")
print("=" * 80)
print()

with open(rule_file, encoding='utf-8') as f:
    data = json.load(f)

changes = []

# Fix default_rule
if 'default_rule' in data and 'conditions' in data['default_rule']:
    cond = data['default_rule']['conditions']
    for field in ['catalyst_loading_molpct', 'base_equiv', 'solvent_volume_mL_per_mmol']:
        if field in cond and isinstance(cond[field], (int, float)):
            old_val = cond[field]
            cond[field] = str(old_val)
            changes.append(f"default_rule.conditions.{field}: {old_val} → \"{old_val}\"")

# Fix base_rules
for idx, rule in enumerate(data.get('base_rules', [])):
    if 'conditions' in rule:
        cond = rule['conditions']
        for field in ['catalyst_loading_molpct', 'base_equiv', 'solvent_volume_mL_per_mmol']:
            if field in cond and isinstance(cond[field], (int, float)):
                old_val = cond[field]
                cond[field] = str(old_val)
                changes.append(f"base_rules[{idx}].conditions.{field}: {old_val} → \"{old_val}\"")

print(f"Found {len(changes)} changes to make:")
print()
for change in changes:
    print(f"  ✓ {change}")

print()
print("=" * 80)
print("PREVIEW OF FIXED FILE")
print("=" * 80)
print()
print("default_rule.conditions (first 5 fields):")
print(json.dumps(
    {k: data['default_rule']['conditions'][k] for k in list(data['default_rule']['conditions'].keys())[:5]},
    indent=2,
    ensure_ascii=False
))

print()
input("Press Enter to save changes, or Ctrl+C to cancel...")

# Save
with open(rule_file, 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=2, ensure_ascii=False)
    f.write('\n')  # Add trailing newline

print()
print("✅ Saved!")
print()
print("Run validation:")
print("  python -m chemtools.schema.builder validate --rules data/rule_db_v2/Suzuki_db.json")
