"""Check motif detection profile for reaction type detection."""
from chemtools.featurizers.detection.matchers import _detect_motif_profile, match_reaction_definition
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

# Test reactants
reactant_smiles = ['Brc1cncnc1', 'OB(O)c1ccoc1']
product_smiles = ['C1(C2=COC=C2)=CN=CN=C1']

print('=== Detected REACTANT motif profile ===')
profile = _detect_motif_profile(reactant_smiles)
for motif_id, info in profile.items():
    print(f"  {motif_id}: count={info['count']}, molecules={info['molecules']}")

print()
print('=== Detected PRODUCT motif profile ===')
product_profile = _detect_motif_profile(product_smiles)
for motif_id, info in product_profile.items():
    print(f"  {motif_id}: count={info['count']}")

# Check which reaction types match
print()
print('=== Reaction type matching ===')
definitions, _ = load_reaction_catalog()

for rx_id in ['Suzuki_miyaura', 'Arylation_Ar_H']:
    if rx_id in definitions:
        defn = definitions[rx_id]
        match = match_reaction_definition(defn, profile, product_profile)
        print(f"{rx_id}:")
        print(f"  reactants required: {defn.reactants}")
        print(f"  products required: {defn.products}")
        if match:
            print(f"  MATCH: slots={match.slot_evidence}, matched={match.matched_slots}/{match.required_slots}")
        else:
            print(f"  NO MATCH")
