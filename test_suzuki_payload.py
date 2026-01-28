"""Test Suzuki detection with the interactive CLI."""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction

# Your reaction
reaction_smiles = "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"

# Featurize with extended options
options = {
    "confirm_coupling_products": True,
    "detailed": True,
    "motif_site_filter": "substituent",
}

payload = featurize_reaction(reaction_smiles, options=options)

print("=" * 80)
print("REACTION PAYLOAD")
print("=" * 80)

# Extract core fields
print(f"\nreaction_smiles: {payload.get('reaction_smiles')}")
print(f"reaction_type: {payload.get('reaction_type')}")
print(f"confidence: {payload.get('confidence')}")
print(f"reaction_key: {payload.get('reaction_key')}")

# Check extended section
extended = payload.get("extended", {})
print(f"\n--- Extended Section ---")
print(f"Has extended: {bool(extended)}")
if extended:
    detection = extended.get("detection", {})
    print(f"\nDetection section: {bool(detection)}")
    if detection:
        matches = detection.get("matches", [])
        print(f"Number of matches: {len(matches)}")
        for i, match in enumerate(matches[:3]):
            print(f"\n  Match {i+1}:")
            print(f"    reaction_type: {match.get('reaction_type')}")
            print(f"    name: {match.get('name')}")
            print(f"    confidence: {match.get('confidence')}")
            print(f"    matched_slots: {match.get('matched_slots')}")
            print(f"    required_slots: {match.get('required_slots')}")
