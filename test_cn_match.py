#!/usr/bin/env python
"""Test if C_N_Coupling matches the Boc-amine reaction."""

from chemtools.featurizers.detection.matchers import _detect_motif_profile, match_reaction_definition
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

rxn = 'O=c1oc2cc(Br)ccc2cc1-c1ccccc1.CC(C)(C)OC(=O)NCCN>>NCCNc1ccc2cc(-c3ccccc3)c(=O)oc2c1'
reactants, products = rxn.split('>>')

# Build profiles
print("1. Building motif profiles...")
reactant_profile = _detect_motif_profile(reactants.split('.'))
product_profile = _detect_motif_profile([products])

print(f"   Reactant motifs: {sorted(reactant_profile.keys())}")
print(f"   Product motifs: {sorted(product_profile.keys())}")

# Load C_N_Coupling
print("\n2. Loading C_N_Coupling definition...")
catalog, _ = load_reaction_catalog()
c_n_def = catalog.get("C_N_Coupling")

print(f"   Electrophile: {c_n_def.reactants['electrophile'].allowed[:3]}...")
print(f"   Nucleophile: {c_n_def.reactants['nucleophile'].allowed[:5]}...")
print(f"   Product: {c_n_def.products['product'].allowed}")

# Check components
print("\n3. Checking match components...")
print(f"   Has HeteroAr-Br? {reactant_profile.get('HeteroAr-Br', {}).get('count', 0) > 0}")
print(f"   Has RCH2-NH2? {reactant_profile.get('RCH2-NH2', {}).get('count', 0) > 0}")
print(f"   Has HeteroAr-NHR? {product_profile.get('HeteroAr-NHR', {}).get('count', 0) > 0}")

# Try matching
print("\n4. Running matcher...")
match = match_reaction_definition(
    definition=c_n_def,
    detected_motifs=reactant_profile,
    detected_products=product_profile,
)

print(f"   Matched? {match is not None}")
if match:
    print(f"   Slots: {match.matched_slots}/{match.required_slots}")
    print(f"   Evidence: {match.slot_evidence}")
else:
    print("   ❌ NO MATCH - but all components are present!")
