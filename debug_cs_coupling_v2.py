from chemtools.featurizers.detection import detect_reaction_types
from chemtools.featurizers.formatters.reaction import extract_reaction_motifs_per_reactant
from chemtools.context import ReactionContext

# Dithiocarbamate + 4-chloroiodobenzene
rxn = "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"

print("=== Dithiocarbamate C-S Coupling Debug v2 ===\n")

# Extract motifs per reactant
ctx = ReactionContext.from_smiles(rxn)
reactant_motifs, product_motifs = extract_reaction_motifs_per_reactant(ctx)

print("Reactant motifs per molecule:")
for i, motifs in enumerate(reactant_motifs):
    print(f"  Reactant {i}: {[m['compound_id'] for m in motifs]}")

print("\nProduct motifs per molecule:")
for i, motifs in enumerate(product_motifs):
    print(f"  Product {i}: {[m['compound_id'] for m in motifs]}")

# Check reaction detection
result = detect_reaction_types(rxn)
print(f"\nDetection result: {len(result.matches)} matches")
for match in result.matches:
    print(f"  {match['reaction_type']} (conf: {match['confidence']:.2f})")
