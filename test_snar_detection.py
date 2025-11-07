"""
Test SNAr detection in the detection engine
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.detection import detect_reaction

print("="*70)
print("Test: SNAr Reaction Detection")
print("="*70)
print()

# User's SNAr reaction
reaction = "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1"
print(f"Reaction: {reaction}")
print()

# Detect with rule-based only (no ML)
print("Detecting reaction family (rule-based)...")
result = detect_reaction(reaction, use_ml=False)

print()
print(f"Family:     {result.get('family')}")
print(f"Confidence: {result.get('confidence')}")
print(f"Method:     {result.get('method')}")
print(f"Status:     {result.get('status')}")
print()

if result.get("details"):
    details = result["details"]
    if "rule_prediction" in details:
        rule_pred = details["rule_prediction"]
        print("Rule prediction:")
        print(f"  Raw family:  {rule_pred.get('raw_family')}")
        print(f"  Confidence:  {rule_pred.get('confidence')}")
        
    if "functional_groups" in details:
        fg = details["functional_groups"]
        print()
        print("Key functional groups detected:")
        relevant = ["aryl_halide", "nucleophile_o", "nucleophile_n", 
                    "catalyst_pd", "catalyst_cu", "catalyst_ni"]
        for key in relevant:
            if key in fg:
                print(f"  {key:25s}: {fg[key]}")

print()
if result.get("family") == "snar":
    print("✅ SUCCESS: SNAr correctly detected!")
else:
    print(f"❌ FAILED: Detected as '{result.get('family')}' instead of 'snar'")
    print()
    print("Full result:")
    import json
    print(json.dumps(result, indent=2))

print("="*70)
