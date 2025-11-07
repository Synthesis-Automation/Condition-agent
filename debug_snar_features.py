"""
Debug feature detection for SNAr
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.featurizers.molecular import featurize

reaction_smiles = "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1"
parts = reaction_smiles.split(">>")
reactants = parts[0].split(".")

print("Electrophile:", reactants[0])
print("Nucleophile:", reactants[1])
print()

features = featurize(reactants[0], reactants[1])

print(f"Total features: {len(features)}")
print()

# Show all True boolean features
print("TRUE boolean features:")
for k, v in sorted(features.items()):
    if isinstance(v, bool) and v:
        print(f"  {k}")

print()

# Show key features the SNAr rule expects
snar_expected = [
    "aromatic_electrophile_present",
    "snar_applicable_electrophile_present",
    "phenol_present",
    "alcohol_present",
    "thiol_present",
    "primary_amine_present",
    "secondary_amine_present",
    "aniline_present"
]

print("SNAr rule file expects these features:")
for feat in snar_expected:
    value = features.get(feat, "NOT FOUND")
    print(f"  {feat:45s}: {value}")
