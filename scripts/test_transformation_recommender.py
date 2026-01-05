import pandas as pd
import os
import sys

# Add current directory to sys.path to import chemtools
sys.path.append(os.getcwd())

from chemtools.HTE.recommender import HTERecommender

def test_recommender():
    db_path = "data/HTE_db/C_N_Coupling_canonical.csv"
    if not os.path.exists(db_path):
        print(f"Error: {db_path} not found. Wait for the conversion script to finish.")
        return

    print(f"Initializing recommender with {db_path}...")
    recommender = HTERecommender(db_path)

    # Test Case 1: Simple Aryl Bromide + Amine
    # Bromobenzene + Aniline
    r1 = "c1ccc(Br)cc1"
    r2 = "c1ccccc1N"
    
    print(f"\n--- Test Case 1: {r1} + {r2} ---")
    result = recommender.recommend(r1, r2)
    
    print(f"Detected Reactant A: {result.reactant_a_type}")
    print(f"Detected Reactant B: {result.reactant_b_type}")
    if result.matched_motifs:
        print(f"Matched Transformation: {result.matched_motifs[0]}")
    else:
        print("No transformation matched.")
    print(f"Total matching experiments: {result.total_matching_experiments}")
    
    if result.recommendations:
        print("\nTop 3 Recommendations:")
        for i, rec in enumerate(result.recommendations[:3]):
            print(f"{i+1}. Catalyst: {rec.catalyst}, Ligand: {rec.ligand}, Base: {rec.base}, Solvent: {rec.solvent}")
            print(f"   Avg Z-score: {rec.avg_z_score:.2f}, Confidence: {rec.confidence_score:.1f}%")
    
    # Test Case 2: Substrate with spectator group
    # 4-Bromophenol + Aniline (Phenol is a spectator in C-N coupling usually)
    r3 = "Oc1ccc(Br)cc1"
    r4 = "c1ccccc1N"
    
    print(f"\n--- Test Case 2: {r3} + {r4} ---")
    result2 = recommender.recommend(r3, r4)
    print(f"Detected Reactant A: {result2.reactant_a_type}")
    if result2.matched_motifs:
        print(f"Matched Transformation: {result2.matched_motifs[0]}")
    else:
        print("No transformation matched.")
    print(f"Total matching experiments: {result2.total_matching_experiments}")
    
    if result2.recommendations:
        print("\nTop 3 Recommendations:")
        for i, rec in enumerate(result2.recommendations[:3]):
            print(f"{i+1}. Catalyst: {rec.catalyst}, Ligand: {rec.ligand}, Base: {rec.base}, Solvent: {rec.solvent}")
            print(f"   Avg Z-score: {rec.avg_z_score:.2f}, Confidence: {rec.confidence_score:.1f}%")

if __name__ == "__main__":
    test_recommender()
