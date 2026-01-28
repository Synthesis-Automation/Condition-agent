"""Test if converter needs update."""

from chemtools.featurizers.formatters.reaction import featurize_reaction
from app.A_convert_to_hte_format import _detect_reaction_type

# User's original Suzuki reaction
smiles = "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"

print("Testing detection paths:\n")

# 1. Converter's detection (used in A_convert_to_hte_format.py)
converter_result = _detect_reaction_type(smiles)
print(f"1. Converter (_detect_reaction_type): {converter_result}")

# 2. Full featurization (with validation)
full_result = featurize_reaction(smiles)
rt = full_result.get("reaction_type", {})
if isinstance(rt, dict):
    rt_name = rt.get("reaction_type", "Unknown")
    rt_conf = rt.get("confidence", 0.0)
else:
    rt_name = str(rt)
    rt_conf = 0.0

print(f"2. Full featurization (with validation): {rt_name} @ {rt_conf}")

# Check if validation was applied
validation = full_result.get("detection", {}).get("validation")
if validation:
    print(f"   ✓ Validation corrected: {validation['original_detection']} → {validation['validated_detection']}")
    print(f"   ✓ Reason: {validation['validation_reason']}")
else:
    print(f"   ℹ No validation correction needed")

print("\n" + "="*60)
if converter_result != rt_name:
    print("⚠️  CONVERTER NEEDS UPDATE")
    print(f"   Converter detects: {converter_result}")
    print(f"   Should detect: {rt_name}")
    print("\nThe converter uses detect_reaction_types() which doesn't include")
    print("validation. It should use the full featurization pipeline to benefit")
    print("from taxonomy-driven validation.")
else:
    print("✅ CONVERTER IS UP TO DATE")
    print("   Detection matches validated result")
