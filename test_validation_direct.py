"""Test validation logic directly"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction

result = featurize_reaction(
    'c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1',
    options={'confirm_coupling_products': True}
)

print("Full reaction bundle:")
print(f"  reaction_type: {result.get('reaction_type')}")
print(f"  confidence: {result.get('confidence')}")

# Check if detection has validation
detection = result.get('detection', {})
print(f"\nDetection payload keys: {list(detection.keys())}")
if 'validation' in detection:
    print(f"Validation: {detection['validation']}")
else:
    print("No validation in detection payload")

# Check aggregates
aggregates = result.get('aggregates')
if aggregates:
    print(f"\nAggregates:")
    print(f"  reacted_motifs: {aggregates.get('reacted_motifs')}")
    print(f"  formed_motifs: {aggregates.get('formed_motifs')}")
