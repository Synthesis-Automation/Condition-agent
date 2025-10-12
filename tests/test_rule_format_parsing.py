"""
Test to verify rule-based format is correctly parsed by LLM synthesis.
Tests the fix for _format_conditions_for_llm to handle chemicals array.
"""

import sys
sys.path.insert(0, '.')

from llmtools.recommendation_llm import _format_conditions_for_llm

print("=" * 70)
print("Testing Rule-Based Format Parsing")
print("=" * 70)

# Rule-based format (from actual output)
rule_conditions = [
    {
        "rank": 1,
        "chemicals": [
            {
                "name": "Tris(dibenzylideneacetone)dipalladium(0)",
                "role": "metal_precursor",
                "abbreviation": "Pd2(dba)3",
                "equivalents": 0.015
            },
            {
                "name": "SPhos",
                "role": "ligand",
                "abbreviation": "SPhos",
                "equivalents": 0.03
            },
            {
                "name": "Potassium carbonate",
                "role": "base",
                "abbreviation": "K2CO3",
                "equivalents": 2.0
            },
            {
                "name": "THF/H2O",
                "role": "solvent",
                "notes": "THF/H2O (4:1)"
            }
        ],
        "conditions": {
            "temperature": [45.0, 60.0],
            "time": [6.0, 10.0]
        }
    }
]

# ML/Protocol flat format (for comparison)
flat_conditions = [
    {
        "catalyst": "Pd(PPh3)4",
        "ligand": "PPh3",
        "solvent": "THF",
        "temperature": "80°C",
        "base": "K2CO3",
        "similarity": 0.95
    }
]

print("\n1. Testing Rule-Based Format (chemicals array):")
print("-" * 70)
rule_text = _format_conditions_for_llm(rule_conditions, "Rule-based")
print(rule_text)

expected_catalyst = "Tris(dibenzylideneacetone)dipalladium(0)"
expected_ligand = "SPhos"
expected_base = "Potassium carbonate"
expected_solvent = "THF/H2O"
expected_temp = "45.0-60.0°C"

if expected_catalyst in rule_text:
    print(f"✅ Catalyst correctly extracted: {expected_catalyst}")
else:
    print(f"❌ Catalyst missing! Expected: {expected_catalyst}")
    sys.exit(1)

if expected_ligand in rule_text:
    print(f"✅ Ligand correctly extracted: {expected_ligand}")
else:
    print(f"❌ Ligand missing! Expected: {expected_ligand}")
    sys.exit(1)

if expected_base in rule_text:
    print(f"✅ Base correctly extracted: {expected_base}")
else:
    print(f"❌ Base missing! Expected: {expected_base}")
    sys.exit(1)

if expected_solvent in rule_text:
    print(f"✅ Solvent correctly extracted: {expected_solvent}")
else:
    print(f"❌ Solvent missing! Expected: {expected_solvent}")
    sys.exit(1)

if expected_temp in rule_text:
    print(f"✅ Temperature correctly extracted: {expected_temp}")
else:
    print(f"❌ Temperature missing! Expected: {expected_temp}")
    sys.exit(1)

print("\n2. Testing Flat Format (backward compatibility):")
print("-" * 70)
flat_text = _format_conditions_for_llm(flat_conditions, "ML-based")
print(flat_text)

if "Pd(PPh3)4" in flat_text and "PPh3" in flat_text and "K2CO3" in flat_text:
    print("✅ Flat format still works correctly")
else:
    print("❌ Flat format broken!")
    sys.exit(1)

print("\n3. Testing Empty Conditions:")
print("-" * 70)
empty_text = _format_conditions_for_llm([], "Empty-test")
print(empty_text)

if "No Empty-test recommendations available" in empty_text:
    print("✅ Empty conditions handled correctly")
else:
    print("❌ Empty conditions not handled!")
    sys.exit(1)

print("\n" + "=" * 70)
print("✅ All format parsing tests passed!")
print("=" * 70)
print("\nRule-based recommendations will now be correctly included in LLM synthesis.")
print("The LLM will see:")
print(f"  - Catalyst: {expected_catalyst}")
print(f"  - Ligand: {expected_ligand}")
print(f"  - Solvent: {expected_solvent}")
print(f"  - Temperature: {expected_temp}")
print(f"  - Base: {expected_base}")
print()
