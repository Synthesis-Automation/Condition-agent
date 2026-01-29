from chemtools.featurizers.detection import detect_reaction_types

# Dithiocarbamate + 4-chloroiodobenzene → aryl dithiocarbamate
rxn = "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"

print("Testing dithiocarbamate C-S coupling:")
print(f"Reaction: {rxn}\n")

result = detect_reaction_types(rxn)

if result.matches:
    top = result.matches[0]
    top_dict = top.to_dict()
    print(f"✓ Detected: {top_dict['reaction_type']} (confidence: {top_dict['confidence']:.2f})")
    print(f"  Nucleophile: {top_dict.get('slots', {}).get('nucleophile')}")
    print(f"  Electrophile: {top_dict.get('slots', {}).get('electrophile')}")
    print(f"  Product: {top_dict.get('slots', {}).get('product')}")
else:
    print("✗ No detection")

print(f"\nAll matches: {len(result.matches)}")
for i, match in enumerate(result.matches[:5]):
    m_dict = match.to_dict()
    print(f"  {i+1}. {m_dict['reaction_type']} (conf: {m_dict['confidence']:.2f})")
