import sys
# Force reimport to clear any caches
if 'chemtools' in sys.modules:
    import importlib
    to_reload = [k for k in sys.modules.keys() if k.startswith('chemtools')]
    for mod in to_reload:
        del sys.modules[mod]

from chemtools.featurizers.detection import detect_reaction_types

# Dithiocarbamate + 4-chloroiodobenzene → aryl dithiocarbamate
rxn = "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"

print("Testing dithiocarbamate C-S coupling (force reload):")
print(f"Reaction: {rxn}\n")

result = detect_reaction_types(rxn)

if result.matches:
    top = result.matches[0]
    print(f"✓ Detected: {top['reaction_type']} (confidence: {top['confidence']:.2f})")
    print(f"  Nucleophile: {top.get('slots', {}).get('nucleophile')}")
    print(f"  Electrophile: {top.get('slots', {}).get('electrophile')}")
    print(f"  Product: {top.get('slots', {}).get('product')}")
else:
    print("✗ No detection")

print(f"\nAll matches: {len(result.matches)}")
