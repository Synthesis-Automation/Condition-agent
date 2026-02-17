"""
Test script to understand HTERecommender workflow and taxonomy system.
"""

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction
from chemtools.featurizers.formatters.reaction import get_crk_options
from chemtools.recommend.recommender import HTERecommender
import json

print("=" * 80)
print("TEST 1: Molecule Featurization (Taxonomy)")
print("=" * 80)

# Test aryl bromide
smiles = "Brc1ccccc1"
result = featurize_molecule(smiles)
print(f"\nSMILES: {smiles}")
print(f"Motifs detected: {[m.get('id') for m in result.get('motifs', [])]}")

# Test boronic acid
smiles = "B(O)(O)c1ccccc1"
result = featurize_molecule(smiles)
print(f"\nSMILES: {smiles}")
print(f"Motifs detected: {[m.get('id') for m in result.get('motifs', [])]}")

print("\n" + "=" * 80)
print("TEST 2: Reaction Featurization")
print("=" * 80)

rxn_smiles = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
result = featurize_reaction(rxn_smiles, options=get_crk_options())
rxn_data = result.get('reaction', {})

print(f"\nReaction SMILES: {rxn_smiles}")
print(f"\nReaction type detected: {rxn_data.get('reaction_type')}")

# Check reactants
reactants = rxn_data.get('reactants', [])
print(f"\nNumber of reactants: {len(reactants)}")
for i, r in enumerate(reactants, 1):
    motifs = [m.get('id') for m in r.get('motifs', [])]
    print(f"  Reactant {i} motifs: {motifs}")

# Check aggregates (reacted, formed, spectator)
aggregates = rxn_data.get('aggregates', {})
print(f"\nReacted motifs: {aggregates.get('reacted_motifs', [])}")
print(f"Formed motifs: {aggregates.get('formed_motifs', [])}")
print(f"Spectator motifs: {aggregates.get('spectator_motifs', [])}")

# Check reaction key
reaction_key = rxn_data.get('reaction_key')
print(f"\nReaction key: {reaction_key}")

print("\n" + "=" * 80)
print("TEST 3: HTERecommender End-to-End")
print("=" * 80)

recommender = HTERecommender('data/HTE_db')

# Test simple Suzuki
reactant_a = "Brc1ccccc1"
reactant_b = "B(O)(O)c1ccccc1"
product = "c1ccc(-c2ccccc2)cc1"

print(f"\nReactant A: {reactant_a}")
print(f"Reactant B: {reactant_b}")
print(f"Product: {product}")

result = recommender.recommend(
    reactant_a_smiles=reactant_a,
    reactant_b_smiles=reactant_b,
    product_smiles=product,
    top_k=5,
    source_group="experiments"
)

print(f"\n--- HTERecommender Internal Detection ---")
print(f"Detected reactant_a_type: {result.reactant_a_type}")
print(f"Detected reactant_b_type: {result.reactant_b_type}")
print(f"Detected product_type: {result.product_type}")
print(f"Predicted reaction type: {result.predicted_reaction_type}")
print(f"Reaction type confidence: {result.reaction_type_confidence:.2f}")
print(f"Reacted motifs: {result.reacted_motifs}")
print(f"Formed motifs: {result.formed_motifs}")
print(f"Spectator motifs: {result.spectator_motifs}")

print(f"\n--- Recommendations ---")
print(f"Found {len(result.recommendations)} recommendations")
for i, rec in enumerate(result.recommendations[:3], 1):
    print(f"\n{i}. Score: {rec.avg_z_score:.2f}")
    print(f"   Catalyst: {rec.catalyst}")
    print(f"   Ligand: {rec.ligand}")
    print(f"   Base: {rec.base}")
    print(f"   Solvent: {rec.solvent}")
    print(f"   Success rate: {rec.success_rate:.1f}%")

print("\n" + "=" * 80)
print("TEST 4: What if we DON'T provide product?")
print("=" * 80)

result_no_product = recommender.recommend(
    reactant_a_smiles=reactant_a,
    reactant_b_smiles=reactant_b,
    product_smiles=None,  # No product
    top_k=5,
    source_group="experiments"
)

print(f"\nDetected reactant_a_type: {result_no_product.reactant_a_type}")
print(f"Detected reactant_b_type: {result_no_product.reactant_b_type}")
print(f"Predicted reaction type: {result_no_product.predicted_reaction_type}")
print(f"Found {len(result_no_product.recommendations)} recommendations")

print("\n" + "=" * 80)
print("TEST 5: Check reaction type mapping")
print("=" * 80)

# Check if Tier 2 names match database
tier2_name = "Suzuki-Miyaura"
print(f"\nTier 2 might say: '{tier2_name}'")
print(f"Database uses: 'Suzuki_miyaura' (note: lowercase 'm')")

# Test with reaction_type_filter
result_filtered = recommender.recommend(
    reactant_a_smiles=reactant_a,
    reactant_b_smiles=reactant_b,
    product_smiles=product,
    top_k=5,
    reaction_type_filter="Suzuki_miyaura",  # Exact database name
    source_group="experiments"
)
print(f"\nWith reaction_type_filter='Suzuki_miyaura': {len(result_filtered.recommendations)} recommendations")

print("\n" + "=" * 80)
print("Summary")
print("=" * 80)
print("""
Key Findings:
1. Taxonomy: featurize_molecule() detects motifs like "Ar-Br", "Ar-B(OR)2"
2. Reaction: featurize_reaction() detects reaction_type and reacted/formed/spectator motifs
3. HTERecommender internally:
   - Calls featurize_molecule() for each reactant
   - Calls featurize_reaction() if product provided
   - Uses detected motifs + reaction_type to match database
4. Database reaction types are like "Suzuki_miyaura" (underscore, specific casing)
5. Product SMILES helps improve matching (reacted/formed motifs more precise)
""")
