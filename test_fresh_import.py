"""Force reimport to test acyl halide detection"""
import importlib
import sys

# Clear chemtools modules from cache
modules_to_clear = [k for k in sys.modules.keys() if k.startswith('chemtools')]
for mod in modules_to_clear:
    del sys.modules[mod]

# Reimport
from chemtools.detection import _DetectionEngine

reaction = "CC(=O)Cl.NCC>>CC(=O)NCC"
e = _DetectionEngine(reaction)
h = e._detect_functional_groups()

print("After fresh import:")
print(f"  acyl_halide: {h.get('acyl_halide')}")
print(f"  nucleophile_n: {h.get('nucleophile_n')}")

result = e._rule_based_detection()
print(f"\nRule detection:")
print(f"  raw_family: {result.get('raw_family')}")
print(f"  family: {result.get('family')}")
