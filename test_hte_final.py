import sys
import os
from pathlib import Path

# Add current directory to path
sys.path.append(os.getcwd())

from chemtools.HTE.recommender import HTERecommender
from examples.sample_reactions import SAMPLE_REACTIONS
from chemtools.smiles import normalize_reaction

def test_hte_recommender_full():
    # Initialize with the V2 database
    recommender = HTERecommender(hte_db_path="data/HTE_db/HTE_V2_Full.csv")
    
    # Test a diverse set of reactions from sample_reactions.py
    # We'll pick a few from different categories
    test_indices = [
        1,   # Suzuki - Simple Ph-Ph
        10,  # Suzuki - Heterocycle
        71,  # Buchwald-Hartwig - Diphenylamine
        72,  # Buchwald-Hartwig - Pyridine ethylamine (Fallback case)
        75,  # C-N Coupling - Indole
        80,  # C-N Coupling - Imidazole
    ]
    
    for idx in test_indices:
        if idx >= len(SAMPLE_REACTIONS):
            continue
            
        rxn_str = SAMPLE_REACTIONS[idx]
        rxn_smiles = rxn_str.split(" (")[0]
        description = rxn_str.split(" (")[1].rstrip(")") if " (" in rxn_str else "No description"
        
        print(f"\n{'='*80}")
        print(f"Testing Reaction {idx}: {description}")
        
        normalized = normalize_reaction(rxn_smiles)
        reactants = normalized.get("reactants") or []
        if len(reactants) < 2:
            print("Skipping: less than 2 reactants")
            continue
            
        smi_a = reactants[0].get("smiles_norm") or reactants[0].get("input")
        smi_b = reactants[1].get("smiles_norm") or reactants[1].get("input")
        
        result = recommender.recommend(smi_a, smi_b)
        
        if result.recommendations:
            print(f"Detected Types: {result.reactant_a_type} + {result.reactant_b_type}")
            if result.is_fallback_match:
                print(f"MATCHED VIA FALLBACK: {result.matched_motifs[0]} + {result.matched_motifs[1]}")
            else:
                print(f"Matched Motifs: {result.matched_motifs[0]} + {result.matched_motifs[1]}")
                
            print(f"Total Experiments in DB: {result.total_matching_experiments}")
            print(f"Predicted Reaction Type: {result.predicted_reaction_type} ({result.reaction_type_confidence*100:.1f}% confidence)")
            
            print("\nTop 3 Recommendations:")
            for i, rec in enumerate(result.recommendations[:3]):
                print(f"{i+1}. Success: {rec.success_rate:.1f}% (Yield: {rec.avg_yield:.1f}%, Z: {rec.avg_z_score:.2f})")
                print(f"   Cat: {rec.catalyst}")
                print(f"   Lig: {rec.ligand}")
                print(f"   Base: {rec.base}")
                print(f"   Sol: {rec.solvent}")
        else:
            print(f"No HTE data found for: {result.reactant_a_type} + {result.reactant_b_type}")

if __name__ == "__main__":
    test_hte_recommender_full()
