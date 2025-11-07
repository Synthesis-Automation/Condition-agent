"""
Debug the featurizer to see why it's not generating the required tokens.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.featurizers.molecular import featurize

# User's SNAr reaction
reaction_smiles = "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1"
reactants = ["Clc1nc(Cl)nc(Cl)n1", "CO"]

print("=" * 70)
print("Debug Featurizer Output")
print("=" * 70)
print()

print(f"Electrophile: {reactants[0]}")
print(f"Nucleophile: {reactants[1]}")
print()

# Featurize
features = featurize(reactants[0], reactants[1])

print("All features:")
for key, value in sorted(features.items()):
    print(f"  {key}: {value}")
print()

# Check for required tokens
required_all = ['aromatic_electrophile_present', 'snar_applicable_electrophile_present']
required_any = ['phenol_present', 'alcohol_present', 'thiol_present', 
                'primary_amine_present', 'secondary_amine_present', 'aniline_present']

print("Required tokens (all):")
for token in required_all:
    present = features.get(token, False)
    status = "✓" if present else "✗"
    print(f"  {status} {token}: {present}")
print()

print("Required tokens (any):")
any_found = False
for token in required_any:
    present = features.get(token, False)
    if present:
        any_found = True
    status = "✓" if present else "✗"
    print(f"  {status} {token}: {present}")
print()

if all(features.get(t) for t in required_all) and any_found:
    print("=" * 70)
    print("✅ All required tokens present!")
    print("=" * 70)
else:
    print("=" * 70)
    print("❌ Missing required tokens")
    print("=" * 70)
