
import pandas as pd
import os
from chemtools.featurizers.reaction_detection import detect_reaction_types
from chemtools.featurizers.unified import featurize_molecule
from rdkit import Chem

def get_motifs(smi):
    """Helper to get unique motif IDs for a SMILES."""
    motifs = set()
    try:
        feat = featurize_molecule(smi)
        if feat and 'molecule' in feat:
            mol_data = feat['molecule']
            if 'motifs' in mol_data:
                for m in mol_data['motifs']:
                    if 'compound_id' in m:
                        motifs.add(m['compound_id'])
    except Exception:
        pass
    return sorted(list(motifs))

def run_comprehensive_detection():
    input_csv = "examples/sample_reactions_2.csv"
    output_csv = "sample_reactions_2_comprehensive_results.csv"
    
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    
    results = []
    print(f"Processing {len(df)} reactions...")
    
    for i, row in df.iterrows():
        # Use 'rxn_smiles_clean' as the column name
        rxn_smiles = row['rxn_smiles_clean']
        
        # 1. Detect Reaction Types
        rxn_res = detect_reaction_types(rxn_smiles)
        detected_rxns = [m.reaction_type for m in rxn_res.matches]
        
        # 2. Split reaction SMILES
        parts = rxn_smiles.split('>')
        reactant_part = parts[0]
        product_part = parts[2] if len(parts) > 2 else ""
        
        reactant_smiles_list = reactant_part.split('.')
        product_smiles_list = product_part.split('.') if product_part else []
        
        # 3. Detect Reactant Motifs (Separately)
        reactant_motifs_list = []
        for smi in reactant_smiles_list:
            motifs = get_motifs(smi)
            reactant_motifs_list.append("|".join(motifs))
        
        # 4. Detect Product Motifs
        product_motifs_list = []
        for smi in product_smiles_list:
            motifs = get_motifs(smi)
            product_motifs_list.append("|".join(motifs))
        
        res_entry = {
            "index": i,
            "reaction_smiles": rxn_smiles,
            "detected_reaction_types": ",".join(detected_rxns),
            "product_motifs": ",".join(product_motifs_list),
            "rxn_match_count": len(detected_rxns),
        }
        
        # Add individual reactant columns
        for idx, motifs_str in enumerate(reactant_motifs_list):
            res_entry[f"reactant_{idx+1}_motifs"] = motifs_str
            
        results.append(res_entry)
        
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(df)}...")

    res_df = pd.DataFrame(results)
    res_df.to_csv(output_csv, index=False)
    print(f"\nDone! Results saved to {output_csv}")
    
    # Print a small summary
    print(f"Total Reactions: {len(res_df)}")
    print(f"Reactions with Type Matches: {len(res_df[res_df['rxn_match_count'] > 0])}")

if __name__ == "__main__":
    run_comprehensive_detection()
