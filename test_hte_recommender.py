"""
Test HTE recommender for the same reaction to compare with precedent-based.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.HTE import HTERecommender

# Test reaction SMILES
reaction_smiles = "c1ccc2c(c1)CCNC2.Nc1ccc(F)c(F)c1>>Nc1ccc(F)c(N2CCc3ccccc3C2)c1"

print(f"Testing HTE Recommender for:")
print(f"Reaction: {reaction_smiles}\n")

# Parse reactants
parts = reaction_smiles.split(">>")
reactants_part = parts[0] if len(parts) > 0 else ""
product = parts[1] if len(parts) > 1 else None

reactants = reactants_part.split(".")
reactant_a = reactants[0] if len(reactants) > 0 else ""
reactant_b = reactants[1] if len(reactants) > 1 else None

print(f"Reactant A: {reactant_a}")
print(f"Reactant B: {reactant_b}")
print(f"Product: {product}\n")

# Initialize HTE recommender
db_path = "data/HTE_db"
print(f"Loading HTE data from {db_path}...")
recommender = HTERecommender(db_path)

print("Running HTE recommendation...\n")
print("=" * 80)

# Run recommendation
result = recommender.recommend(
    reactant_a_smiles=reactant_a,
    reactant_b_smiles=reactant_b,
    product_smiles=product,
    top_k=10,
    min_experiments=1,
    reaction_type_filter=None,
    catalyst_filter=None,
    source_group=None,
    use_aryl_steric_electronic_weighting=False,
)

print(f"Reactant A Type: {result.reactant_a_type}")
print(f"Reactant A Category: {result.reactant_a_category}")
print(f"Reactant B Type: {result.reactant_b_type}")
print(f"Reactant B Category: {result.reactant_b_category}")
if result.product_type:
    print(f"Product Type: {result.product_type}")
if result.reacted_motifs:
    print(f"Reacted Motifs: {result.reacted_motifs}")
if result.formed_motifs:
    print(f"Formed Motifs: {result.formed_motifs}")
print()

# Show recommendations by source
for source_type, recommendations in result.recommendations_by_source.items():
    print(f"\n{source_type.upper()} Recommendations: {len(recommendations)} found")
    print("-" * 80)
    
    for i, rec in enumerate(recommendations[:10], 1):
        print(f"\nRank {i}:")
        print(f"  Reaction Type: {rec.reaction_type}")
        if hasattr(rec, 'reaction_id') and rec.reaction_id:
            print(f"  Reaction ID: {rec.reaction_id}")
        if hasattr(rec, 'match_score'):
            print(f"  Match Score: {rec.match_score:.3f}")
        if hasattr(rec, 'avg_z_score'):
            print(f"  Avg Z-Score: {rec.avg_z_score:.3f}")
        print(f"  Catalyst: {rec.catalyst or 'None'}")
        print(f"  Ligand: {rec.ligand or 'None'}")
        print(f"  Base: {rec.base or 'None'}")
        print(f"  Solvent: {rec.solvent or 'None'}")
        print(f"  Additive: {rec.additive or 'None'}")
        if hasattr(rec, 'avg_yield') and rec.avg_yield is not None:
            print(f"  Avg Yield: {rec.avg_yield:.1f}%")
        if hasattr(rec, 'success_rate') and rec.success_rate is not None:
            print(f"  Success Rate: {rec.success_rate:.1f}%")
        if hasattr(rec, 'num_experiments'):
            print(f"  Experiments: {rec.num_experiments}")

print("\n" + "=" * 80)
print("Done!")
