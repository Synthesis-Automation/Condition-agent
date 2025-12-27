
import sys
from pathlib import Path

# Add the current directory to sys.path to import chemtools
sys.path.append(str(Path.cwd()))

from chemtools.featurizers.structural import featurize_molecule
from examples.sample_compounds import ALL_SAMPLE_COMPOUNDS
from rdkit import Chem

def run_test():
    # Manual test cases for new types
    MANUAL_TEST_CASES = [
        {"name": "Isopropanol", "smiles": "CC(C)O"},
        {"name": "tert-Butanol", "smiles": "CC(C)(C)O"},
        {"name": "Acetamide", "smiles": "CC(=O)N"},
        {"name": "Methanesulfonamide", "smiles": "CS(=O)(=O)N"},
        {"name": "Pyrrole", "smiles": "c1cc[nH]c1"},
        {"name": "Diphenylamine", "smiles": "c1ccc(Nc2ccccc2)cc1"},
        {"name": "2-Piperidone", "smiles": "O=C1CCCCN1"},
        {"name": "N-Chloropyrrole", "smiles": "Cln1cccc1"},
        {"name": "4-Bromopyridine", "smiles": "Brc1ccncc1"},
        {"name": "Phenylphosphonium", "smiles": "c1ccc([P+](c2ccccc2)(c3ccccc3)c4ccccc4)cc1"},
        {"name": "Phenylboronate", "smiles": "c1ccc(B(OC)OC)cc1"},
    ]

    print(f"Running analysis on {len(ALL_SAMPLE_COMPOUNDS)} sample compounds + {len(MANUAL_TEST_CASES)} manual cases...\n")
    
    motif_counts = {}
    
    # We'll look for specific interesting compounds to show detailed results
    interesting_names = [c["name"] for c in MANUAL_TEST_CASES]
    
    detailed_results = []

    for comp in ALL_SAMPLE_COMPOUNDS + MANUAL_TEST_CASES:
        smiles = comp["smiles"]
        name = comp.get("name", "Unknown")
        
        try:
            # Try with AddHs to see if H-motifs show up
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol_with_hs = Chem.AddHs(mol)
                # Note: featurize_molecule takes SMILES, so it will parse it again.
                # To test AddHs, we'd need to modify featurize_molecule or call detect_motifs directly.
                # For now, let's just see what the standard featurize_molecule does.
                
                analysis = featurize_molecule(smiles)
                motifs = [m["compound_id"] for m in analysis.get("motifs", [])]
                
                for m in motifs:
                    motif_counts[m] = motif_counts.get(m, 0) + 1
                
                if name in interesting_names:
                    detailed_results.append({
                        "name": name,
                        "smiles": smiles,
                        "detected": motifs,
                        "analysis": analysis.get("analyses", [])
                    })
        except Exception as e:
            print(f"Error analyzing {name} ({smiles}): {e}")

    print("--- Motif Detection Summary ---")
    # Sort by count descending
    sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)
    for motif, count in sorted_motifs:
        print(f"{motif:20}: {count}")
    
    print("\n--- Detailed Results for Selected Compounds ---")
    for res in detailed_results:
        print(f"Name: {res['name']}")
        print(f"SMILES: {res['smiles']}")
        print(f"Detected Motifs: {res['detected']}")
        if res['analysis']:
            for a in res['analysis']:
                cid = a['compound_id']
                steric = a.get('steric', {}).get('score_0_10', 'N/A')
                elec = a.get('electronic', {}).get('score_0_10', 'N/A')
                print(f"  - {cid}: Steric={steric}, Electronic={elec}")
        print("-" * 40)

if __name__ == "__main__":
    run_test()
