"""Find files with actual mixed types that need fixing."""

import json
from pathlib import Path
from collections import defaultdict

rule_dir = Path('data/rule_db_v2')

issues_by_file = defaultdict(list)
numeric_fields = {
    "catalyst_loading_molpct",
    "ligand_loading_molpct",
    "pd_loading_molpct",
    "cu_loading_molpct",
    "base_equiv",
    "temperature_C",
    "time_h",
    "pressure_bar",
    "concentration_M",
    "substrate_concentration_mM",
    "solvent_volume_mL_per_mmol",
}

print("=" * 80)
print("FINDING ACTUAL MIXED TYPE ISSUES")
print("=" * 80)
print()

for rule_file in sorted(rule_dir.glob('*.json')):
    with open(rule_file, encoding='utf-8') as f:
        data = json.load(f)
    
    # Check default_rule
    default = data.get('default_rule', {}).get('conditions', {})
    for field, value in default.items():
        if field in numeric_fields and not isinstance(value, str):
            issues_by_file[rule_file.name].append({
                'location': 'default_rule.conditions',
                'field': field,
                'type': type(value).__name__,
                'value': value
            })
    
    # Check base_rules
    for idx, br in enumerate(data.get('base_rules', [])):
        conditions = br.get('conditions', {})
        for field, value in conditions.items():
            if field in numeric_fields and not isinstance(value, str):
                issues_by_file[rule_file.name].append({
                    'location': f'base_rules[{idx}].conditions',
                    'field': field,
                    'type': type(value).__name__,
                    'value': value
                })

if not issues_by_file:
    print("✅ NO MIXED TYPE ISSUES FOUND!")
    print()
    print("Your rule files are already using string types consistently.")
    print("The standardization focus should be on:")
    print("  1. Field naming consistency")
    print("  2. Documentation")
    print("  3. Validation rules")
else:
    total_issues = sum(len(issues) for issues in issues_by_file.values())
    print(f"Found {total_issues} mixed type issues across {len(issues_by_file)} files:")
    print()
    
    for filename, issues in sorted(issues_by_file.items()):
        print(f"📁 {filename}: {len(issues)} issues")
        for issue in issues:
            print(f"   {issue['location']}.{issue['field']}")
            print(f"   Current: {issue['type']} = {issue['value']}")
            print(f"   Fix: Change to string \"{issue['value']}\"")
            print()

print()
print("=" * 80)
print("RECOMMENDED ACTIONS")
print("=" * 80)
print()

if not issues_by_file:
    print("Since types are already consistent, focus on:")
    print()
    print("1. Create SCHEMA_CONDITIONS.md (documentation)")
    print("   - Define standard field names")
    print("   - Document value formats")
    print("   - List family-specific fields")
    print()
    print("2. Enhance validator.py")
    print("   - Add field name checks")
    print("   - Add range format validation")
    print("   - Add required field checks")
    print()
    print("3. Proceed directly to Addition Sequence implementation")
    print("   - Your data is clean enough to use now")
    print("   - Can refine as you go")
    print()
    print("⏱️ Time saved: ~6 hours (no file modifications needed)")
else:
    print("1. Run auto-fix script to convert numbers to strings")
    print("2. Validate changes")
    print("3. Proceed with addition sequence")
