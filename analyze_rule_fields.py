"""Analyze condition fields across all rule files to identify standardization needs."""

import json
from pathlib import Path
from collections import Counter, defaultdict

rule_dir = Path('data/rule_db_v2')

# Collect all field names and their variations
all_fields = Counter()
field_examples = defaultdict(set)
field_types = defaultdict(set)

for rule_file in sorted(rule_dir.glob('*.json')):
    with open(rule_file, encoding='utf-8') as f:
        data = json.load(f)
    
    # Check default_rule
    default = data.get('default_rule', {}).get('conditions', {})
    for key, value in default.items():
        all_fields[key] += 1
        field_examples[key].add(str(value)[:80])
        field_types[key].add(type(value).__name__)
    
    # Check base_rules
    for br in data.get('base_rules', []):
        conditions = br.get('conditions', {})
        for key, value in conditions.items():
            all_fields[key] += 1
            field_examples[key].add(str(value)[:80])
            field_types[key].add(type(value).__name__)

print("=" * 80)
print("CONDITION FIELDS FOUND ACROSS ALL RULE FILES")
print("=" * 80)
print(f"\nTotal unique fields: {len(all_fields)}")
print(f"Total rule files analyzed: {len(list(rule_dir.glob('*.json')))}")

print("\n" + "=" * 80)
print("FIELD USAGE FREQUENCY")
print("=" * 80)
for field, count in sorted(all_fields.items(), key=lambda x: -x[1]):
    types_str = ", ".join(sorted(field_types[field]))
    print(f"\n{field}: {count} occurrences (types: {types_str})")
    examples = list(field_examples[field])[:3]
    for ex in examples:
        print(f"  Example: {ex}")

print("\n" + "=" * 80)
print("STANDARDIZATION ISSUES DETECTED")
print("=" * 80)

issues = []

# Check for inconsistent naming
naming_variants = [
    ("catalyst_loading_molpct", "catalyst_loading_mol_pct", "catalyst_mol_pct"),
    ("base_equiv", "base_equivalents"),
    ("ligand_loading_molpct", "ligand_loading_mol_pct", "ligand_mol_pct"),
    ("temperature_C", "temp_C", "temperature"),
    ("time_h", "time_hours", "reaction_time_h"),
]

for variants in naming_variants:
    found = [v for v in variants if v in all_fields]
    if len(found) > 1:
        issues.append(f"⚠️  Inconsistent naming: {', '.join(found)}")
        print(f"\n⚠️  Inconsistent naming:")
        for v in found:
            print(f"    '{v}': {all_fields[v]} occurrences")

# Check for fields that should always be present
required_fields = ["catalyst", "solvent", "temperature_C"]
for field in required_fields:
    if all_fields[field] < 45:  # Assuming ~5 rules per file * 9 files = 45
        issues.append(f"⚠️  '{field}' not consistently present ({all_fields[field]} occurrences)")
        print(f"\n⚠️  '{field}' not consistently present ({all_fields[field]} occurrences)")

# Check for mixed types
for field, types in field_types.items():
    if len(types) > 1:
        issues.append(f"⚠️  '{field}' has mixed types: {', '.join(types)}")
        print(f"\n⚠️  '{field}' has mixed types: {', '.join(types)}")

if not issues:
    print("\n✅ No major standardization issues detected!")
else:
    print(f"\n\nTotal issues found: {len(issues)}")

print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)
print("""
1. Standardize field names (choose one variant):
   - catalyst_loading_molpct (recommended - most common)
   - ligand_loading_molpct (recommended - most common)
   - base_equiv (recommended - most common)
   - temperature_C (recommended - most common)
   - time_h (recommended - most common)

2. Define required vs optional fields in schema

3. Standardize value formats:
   - Ranges: Always use "min-max" format (e.g., "0.5-2.0")
   - Options: Always use " or " separator (e.g., "THF or toluene")
   - Lists: Always use arrays for multiple items

4. Add field-level documentation/schema
""")
