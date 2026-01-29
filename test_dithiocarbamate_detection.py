from chemtools.featurizers.detection import detect_reaction_types

# Dithiocarbamate + 4-chloroiodobenzene → aryl dithiocarbamate
rxn = "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"

print("Testing dithiocarbamate C-S coupling:")
print(f"Reaction: {rxn}\n")

result = detect_reaction_types(rxn)

if result.matches:
    top = result.matches[0]
    print(f"✓ Detected: {top['reaction_type']} (confidence: {top['confidence']:.2f})")
else:
    print("✗ No detection")

print(f"\nAll matches ({len(result.matches)}):")
for match in result.matches:
    print(f"  {match['reaction_type']} (conf: {match['confidence']:.2f})")
    print(f"    Nucleophile: {match.get('slots', {}).get('nucleophile')}")
    print(f"    Electrophile: {match.get('slots', {}).get('electrophile')}")
    print(f"    Product: {match.get('slots', {}).get('product')}")
