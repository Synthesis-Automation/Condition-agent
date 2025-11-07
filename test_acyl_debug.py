"""Debug acyl halide detection"""
from chemtools.detection import _DetectionEngine

reaction = "CC(=O)Cl.NCC>>CC(=O)NCC"
e = _DetectionEngine(reaction)
h = e._detect_functional_groups()

print("Functional groups:")
print(f"  acyl_halide: {h.get('acyl_halide')}")
print(f"  nucleophile_n: {h.get('nucleophile_n')}")
print(f"  nucleophile_o: {h.get('nucleophile_o')}")
print(f"  carbonyl: {h.get('carbonyl')}")

# Check the condition
if h.get("acyl_halide") and (h.get("nucleophile_n") or h.get("nucleophile_o")):
    print("\n✅ Condition matches! Should set fam='amide_formation'")
else:
    print("\n❌ Condition doesn't match")

result = e._rule_based_detection()
print(f"\nRule detection result:")
print(f"  raw_family: {result.get('raw_family')}")
print(f"  family: {result.get('family')}")
print(f"  confidence: {result.get('confidence')}")
