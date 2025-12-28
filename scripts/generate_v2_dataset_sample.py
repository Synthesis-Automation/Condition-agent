import json
from chemtools.featurizers.structural import featurize_reaction
from pathlib import Path

def process_v2(input_path, output_path, limit=10):
    with open(input_path, 'r', encoding='utf-8') as f_in, \
         open(output_path, 'w', encoding='utf-8') as f_out:
        
        for i, line in enumerate(f_in):
            if i >= limit:
                break
            
            data = json.loads(line)
            reaction_smiles = data.get("precomputed", {}).get("reaction_smiles")
            
            if not reaction_smiles:
                continue
            
            # Run new analysis
            try:
                analysis = featurize_reaction(reaction_smiles)
                
                # Create a simplified V2 structure
                v2_entry = {
                    "reaction_id": data.get("reaction_id"),
                    "reaction_type": data.get("reaction_type"),
                    "condition_core": data.get("condition_core"),
                    "catalytic_system": data.get("catalytic_system"),
                    "reagents": data.get("reagents"),
                    "solvents": data.get("solvents"),
                    "conditions": data.get("conditions"),
                    "reference": data.get("reference"),
                    "reaction_smiles": reaction_smiles,
                    "analysis_v2": {
                        "detection": analysis.get("detection"),
                        "motifs": [
                            {
                                "id": m["compound_id"],
                                "smiles": m["smiles"],
                                "atoms": m["bond"]
                            }
                            for m in analysis.get("motifs", [])
                        ],
                        "features": {
                            "steric": analysis.get("steric"),
                            "electronics": analysis.get("electronics"),
                            "nearby": analysis.get("nearby")
                        }
                    }
                }
                f_out.write(json.dumps(v2_entry) + "\n")
                print(f"Processed {i+1}/{limit}")
            except Exception as e:
                print(f"Error processing {i+1}: {e}")

if __name__ == "__main__":
    input_file = r"c:\Git-softwares\Condition-agent\data\reaction_dataset\C_N_Coupling.jsonl"
    output_file = r"c:\Git-softwares\Condition-agent\data\reaction_dataset\C_N_Coupling_V2_Sample.jsonl"
    process_v2(input_file, output_file, limit=5)
