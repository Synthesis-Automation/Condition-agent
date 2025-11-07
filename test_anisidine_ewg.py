"""
Test to see why p-anisidine is being identified as having an EWG.
"""

from chemtools.featurizers.molecular import featurize

# p-Anisidine: para-methoxyaniline (has EDG -OMe, not EWG)
p_anisidine = "Nc1ccc(OC)cc1"

# Test as nucleophile with benzaldehyde
benzaldehyde = "O=Cc1ccccc1"

print("=" * 70)
print("Testing p-Anisidine EWG Detection")
print("=" * 70)
print()
print(f"p-Anisidine SMILES: {p_anisidine}")
print(f"Electrophile: {benzaldehyde}")
print()

# Featurize the pair
features = featurize(benzaldehyde, p_anisidine)

print("Key features:")
print(f"  elec_class: {features.get('elec_class', 'NOT FOUND')}")
print(f"  para_EWG: {features.get('para_EWG', 'NOT FOUND')}")
print(f"  LG: {features.get('LG', 'NOT FOUND')}")
print(f"  nuc_class: {features.get('nuc_class', 'NOT FOUND')}")
print()

# Check what groups are detected
print("Functional groups detected:")
for key, value in sorted(features.items()):
    if "_present" in key and value:
        print(f"  {key}: {value}")
print()

# The issue: methoxy (-OMe) is an EDG by resonance, NOT an EWG
print("Analysis:")
print("  -OMe (methoxy) is an electron-DONATING group (EDG)")
print("  It donates electrons through resonance (+M effect)")
print("  It should NOT be flagged as para_EWG")
print()
print("=" * 70)
