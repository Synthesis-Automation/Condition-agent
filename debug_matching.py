#!/usr/bin/env python
"""Debug: Step through the matching process for C-O coupling."""

from chemtools.featurizers.detection.matchers import _detect_motif_profile, match_reaction_definition
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

# Phenol + iodobenzonitrile → aryl ether
rxn_smiles = 'Cc1cc(C)cc(O)c1.N#Cc1ccc(I)cc1>>Cc1cc(C)cc(Oc2ccc(C#N)cc2)c1'
reactants, products = rxn_smiles.split('>>')

# Build motif profiles
print("1. Building motif profiles...")
reactant_profile = _detect_motif_profile(reactants.split('.'))
product_profile = _detect_motif_profile([products])

print(f"   Reactant motifs: {sorted(reactant_profile.keys())}")
print(f"   Product motifs: {sorted(product_profile.keys())}")

# Check C_O_Coupling definition
print("\n2. Loading C_O_Coupling definition...")
catalog, _ = load_reaction_catalog()
c_o_def = catalog.get("C_O_Coupling")

print(f"   Electrophile required: {c_o_def.reactants['electrophile'].allowed[:3]}...")
print(f"   Nucleophile required: {c_o_def.reactants['nucleophile'].allowed[:3]}...")
print(f"   Product required: {c_o_def.products['product'].allowed}")

# Check matches
print("\n3. Checking matches...")
print(f"   Has Ar-I (electrophile)? {reactant_profile.get('Ar-I', {}).get('count', 0) > 0}")
print(f"   Has Ar-OH (nucleophile)? {reactant_profile.get('Ar-OH', {}).get('count', 0) > 0}")
print(f"   Has Ar-OR (product)? {product_profile.get('Ar-OR', {}).get('count', 0) > 0}")

# Try matching
print("\n4. Running matcher...")
match = match_reaction_definition(
    definition=c_o_def,
    detected_motifs=reactant_profile,
    detected_products=product_profile,
)

print(f"   Matched? {match is not None}")
if match:
    print(f"   Reaction type: {match.reaction_type}")
    print(f"   Evidence: {match.slot_evidence}")
    print(f"   Slots: {match.matched_slots}/{match.required_slots}")
