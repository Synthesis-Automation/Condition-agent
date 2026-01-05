import pandas as pd
from chemtools.HTE.recommender import HTERecommender
from pathlib import Path
import os

def generate_report():
    # Paths
    db_path = "data/HTE_db/C_N_Coupling_canonical.csv"
    samples_path = "examples/sample_reactions.csv"
    output_path = "results/hte_cn_report.csv"
    
    os.makedirs("results", exist_ok=True)
    
    if not os.path.exists(db_path):
        print(f"Error: {db_path} not found")
        return
    
    # Initialize recommender
    print(f"Loading HTE database from {db_path}...")
    recommender = HTERecommender(hte_db_path=db_path)
    
    # Read samples
    print(f"Reading samples from {samples_path}...")
    samples_df = pd.read_csv(samples_path)
    
    # Filter for Buchwald-Hartwig reactions
    cn_samples = samples_df[samples_df['reaction_type'].str.contains('Buchwald-Hartwig', na=False, case=False)]
    print(f"Found {len(cn_samples)} C-N formation reactions.")
    
    report_data = []
    
    for idx, row in cn_samples.iterrows():
        reactants = str(row['reactants_smiles']).split('.')
        if len(reactants) < 2:
            continue
            
        r_a = reactants[0]
        r_b = reactants[1]
        
        print(f"Processing: {row['example']} ({r_a} + {r_b})")
        
        try:
            result = recommender.recommend(r_a, r_b, top_k=1)
            
            if result.recommendations:
                rec = result.recommendations[0]
                report_data.append({
                    "Example": row['example'],
                    "Reactant_A": r_a,
                    "Reactant_B": r_b,
                    "Detected_A": result.reactant_a_type,
                    "Detected_B": result.reactant_b_type,
                    "Predicted_Type": result.predicted_reaction_type,
                    "Confidence": result.reaction_type_confidence,
                    "Matching_Exps": result.total_matching_experiments,
                    "Top_Catalyst": rec.catalyst,
                    "Top_Ligand": rec.ligand,
                    "Top_Base": rec.base,
                    "Top_Solvent": rec.solvent,
                    "Avg_Z_Score": rec.avg_z_score,
                    "Success_Rate": rec.success_rate
                })
            else:
                report_data.append({
                    "Example": row['example'],
                    "Reactant_A": r_a,
                    "Reactant_B": r_b,
                    "Detected_A": result.reactant_a_type,
                    "Detected_B": result.reactant_b_type,
                    "Predicted_Type": "No Match",
                    "Confidence": 0,
                    "Matching_Exps": 0,
                    "Top_Catalyst": "N/A",
                    "Top_Ligand": "N/A",
                    "Top_Base": "N/A",
                    "Top_Solvent": "N/A",
                    "Avg_Z_Score": 0,
                    "Success_Rate": 0
                })
        except Exception as e:
            print(f"Error processing {row['example']}: {e}")
            
    # Save report
    report_df = pd.DataFrame(report_data)
    report_df.to_csv(output_path, index=False)
    print(f"Report saved to {output_path}")

if __name__ == "__main__":
    generate_report()
