"""Debug script with more details."""

from chemtools.detection import _DetectionEngine

rcm_reaction = "C=CCNC(=O)C=C>>C1=CCNC(=O)C1"

print("Testing RCM reaction detection with detailed output:")
print(f"Reaction: {rcm_reaction}\n")

detector = _DetectionEngine(rcm_reaction)

print("Functional groups detected:")
for key, value in detector.functional_groups.items():
    if value:
        print(f"  ✓ {key}: {value}")

print(f"\nReactants: {detector.reactants}")
print(f"Catalysts: {detector.catalysts}")
print(f"Agents: {detector.agents}")

result = detector.detect(use_ml=False)
print(f"\nFinal detection:")
print(f"  Family: {result.get('family')}")
print(f"  Confidence: {result.get('confidence')}")
