"""
Test Reductive Amination rule integration
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.rule.engine import RuleEngine
from chemtools.featurizers.molecular import featurize

print("="*70)
print("Reductive Amination Integration Test")
print("="*70)
print()

# Example reductive amination: aldehyde + amine
electrophile = "CC(=O)C"  # Acetone (ketone)
nucleophile = "CCN"       # Ethylamine (primary amine)

print(f"Electrophile (ketone): {electrophile}")
print(f"Nucleophile (amine):   {nucleophile}")
print()

# Featurize
print("Step 1: Featurizing...")
features = featurize(electrophile, nucleophile)
print(f"✓ Generated {len(features)} features")
print()

# Show key reductive amination features
print("Key reductive amination features:")
ra_keys = [
    "carbonyl_present",
    "ketone_present",
    "aldehyde_present",
    "amine_nucleophile_present",
    "primary_amine_present",
    "amine_primary_present"
]
for key in ra_keys:
    value = features.get(key, "NOT FOUND")
    status = "✓" if value else "✗"
    print(f"  {status} {key:35s}: {value}")
print()

# Load reductive amination rules
print("Step 2: Loading reductive amination rules...")
try:
    engine = RuleEngine.from_file("data/rule_db/reductive_amination_db.json")
    print(f"✓ Loaded reductive_amination_db.json")
    print(f"  Name: {engine.database.metadata.get('name')}")
    print(f"  Rules: {len(engine.database.base_rules)}")
    print()
except Exception as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Check if rules apply
print("Step 3: Checking if reductive amination rules apply...")
applies = engine.database.check_applies(features)
print(f"  Result: {applies}")
if not applies:
    print("  ✗ Rules do not apply")
    print()
    print("  Required by applies_if:")
    for key, value in engine.database.applies_if.items():
        print(f"    {key}: {value}")
    print()
    print("  Present features:")
    for key, value in sorted(features.items()):
        if isinstance(value, bool) and value:
            print(f"    ✓ {key}")
    sys.exit(1)
print(f"  ✓ Reductive amination rules apply!")
print()

# Find matching rule
print("Step 4: Finding matching rule...")
matched_rule = engine.database.find_matching_rule(features)
if not matched_rule:
    print("  Using default rule")
    matched_rule = engine.database.find_matching_rule(features, use_default=True)

if matched_rule:
    print(f"  ✓ Matched: {matched_rule.name}")
    print()
    print("  Key conditions:")
    key_conds = ["reducing_agent", "solvent", "temperature_C", "time_h", "acid_or_buffer"]
    for key in key_conds:
        value = matched_rule.conditions.get(key, "N/A")
        if value != "N/A":
            print(f"    {key:25s}: {value}")
print()

print("="*70)
print("✅ SUCCESS: Reductive Amination Integration Working!")
print("="*70)
