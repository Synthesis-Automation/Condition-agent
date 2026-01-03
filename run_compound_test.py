import pandas as pd
import os
from chemtools.featurizers.unified import featurize_molecule

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

def run_compound_classification():
    input_csv = "examples/sample_compounds.csv"
    output_csv = "sample_compounds_results.csv"
    
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    
    results = []
    print(f"Processing {len(df)} compounds...")
    
    for i, row in df.iterrows():
        smi = row['smiles']
        name = row['name']
        
        motifs = get_motifs(smi)
        
        res_entry = {
            "index": i,
            "name": name,
            "smiles": smi,
            "detected_motifs": "|".join(motifs),
            "motif_count": len(motifs)
        }
        results.append(res_entry)
        
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(df)}...")

    res_df = pd.DataFrame(results)
    res_df.to_csv(output_csv, index=False)
    print(f"\nDone! Results saved to {output_csv}")

if __name__ == "__main__":
    run_compound_classification()
