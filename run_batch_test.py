import csv
import json
from pathlib import Path
from chemtools.featurizers.unified import featurize_reaction

def process_csv(input_path: str, output_path: str):
    input_file = Path(input_path)
    output_file = Path(output_path)
    
    results = []
    with input_file.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles = row.get("reaction_smiles")
            expected = row.get("reaction_type")
            notation = row.get("notation")
            
            if not smiles:
                continue
                
            try:
                payload = featurize_reaction(smiles)
                matches = payload['reaction']['detection']['matches']
                detected_ids = [m.get('reaction_type') for m in matches]
                detected_names = [m.get('name') for m in matches]
                
                results.append({
                    "notation": notation,
                    "expected": expected,
                    "detected_ids": "|".join(detected_ids),
                    "detected_names": "|".join(detected_names),
                    "smiles": smiles
                })
            except Exception as e:
                results.append({
                    "notation": notation,
                    "expected": expected,
                    "detected_ids": "ERROR",
                    "detected_names": str(e),
                    "smiles": smiles
                })

    fieldnames = ["notation", "expected", "detected_ids", "detected_names", "smiles"]
    with output_file.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Processed {len(results)} reactions. Results saved to {output_path}")

if __name__ == "__main__":
    process_csv(
        "c:/Git-softwares/Condition-agent/examples/sample_reactions.csv",
        "c:/Git-softwares/Condition-agent/sample_reactions_results.csv"
    )
