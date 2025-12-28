import json
from chemtools.featurizers.structural import featurize_reaction
from pathlib import Path
import sys

# Add current dir to path to ensure chemtools is found
sys.path.append(".")

def extract_recommendation_features(analysis):
    """
    Extract a flattened dictionary of features most useful for KNN/Recommendation.
    """
    features = {}
    
    # 1. Reaction Type / Family
    detection = analysis.get("detection", {})
    matches = detection.get("matches", [])
    if matches:
        # Use the first/best match
        best_match = matches[0]
        features["reaction_type"] = best_match.get("reaction_type")
        features["reaction_category"] = best_match.get("category")
    
    # 2. Steric Summary
    steric = analysis.get("steric", {})
    # Get max steric score for aryl and alkyl
    aryl_sterics = [s.get("result", {}).get("score_0_10", 0) for s in steric.get("aryl", [])]
    alkyl_sterics = [s.get("result", {}).get("score_0_10", 0) for s in steric.get("alkyl", [])]
    features["max_aryl_steric"] = max(aryl_sterics) if aryl_sterics else 0
    features["max_alkyl_steric"] = max(alkyl_sterics) if alkyl_sterics else 0
    
    # 3. Electronic Summary
    electronics = analysis.get("electronics", {})
    aryl_elec = [e.get("result", {}).get("score_0_10", 5.0) for e in electronics.get("aryl", [])]
    # Average electronic score (5.0 is neutral)
    features["avg_aryl_electronic"] = sum(aryl_elec) / len(aryl_elec) if aryl_elec else 5.0
    
    # 4. Nearby Groups / Chelators
    nearby = analysis.get("nearby", [])
    has_chelator = False
    for n in nearby:
        for res in n.get("result", []):
            tags = res.get("tags", [])
            if "bidentate_chelator" in tags or "chelator" in tags:
                has_chelator = True
                break
    features["has_chelator"] = has_chelator
    
    # 5. Motif Counts (Simplified)
    motifs = analysis.get("motifs", [])
    motif_ids = [m["compound_id"] for m in motifs]
    for mid in set(motif_ids):
        features[f"has_{mid}"] = True
        
    return features

def process_recommendation_dataset(input_path, output_path, limit=-1):
    print(f"Starting migration: {input_path} -> {output_path}")
    
    # Count total lines for progress
    total_lines = 0
    with open(input_path, "r", encoding="utf-8") as f:
        for _ in f: total_lines += 1
    
    if limit > 0:
        total_lines = min(total_lines, limit)
        
    print(f"Total reactions to process: {total_lines}")

    with open(input_path, "r", encoding="utf-8") as f_in, \
         open(output_path, "w", encoding="utf-8") as f_out:
        
        for i, line in enumerate(f_in):
            if limit > 0 and i >= limit:
                break
            
            if i % 100 == 0:
                print(f"Progress: {i}/{total_lines} ({i/total_lines*100:.1f}%)")
            
            try:
                data = json.loads(line)
                reaction_smiles = data.get("precomputed", {}).get("reaction_smiles")
                if not reaction_smiles:
                    continue
                
                # 1. Structural Analysis
                analysis = featurize_reaction(reaction_smiles)
                
                # 2. Recommendation Enrichment
                rec_features = extract_recommendation_features(analysis)
                
                # 3. Condition Normalization
                conditions = {
                    "catalyst_ligand": data.get("condition_core"),
                    "reagents": [r.get("name") for r in data.get("reagents", [])],
                    "solvents": [s.get("name") for s in data.get("solvents", [])],
                    "temp_c": data.get("conditions", {}).get("temperature_c"),
                    "yield": data.get("conditions", {}).get("yield_pct")
                }
                
                # 4. Final V2 Recommendation Entry
                v2_entry = {
                    "reaction_id": data.get("reaction_id"),
                    "reaction_smiles": reaction_smiles,
                    "conditions": conditions,
                    "recommendation_payload": {
                        "family": rec_features.get("reaction_type"),
                        "features": rec_features,
                        "structural_analysis": {
                            "detection": analysis.get("detection"),
                            "steric": analysis.get("steric"),
                            "electronics": analysis.get("electronics"),
                            "nearby": analysis.get("nearby")
                        }
                    },
                    "reference": data.get("reference")
                }
                
                f_out.write(json.dumps(v2_entry) + "\n")
                
            except Exception as e:
                print(f"Error processing line {i}: {e}")
                continue

    print(f"Migration complete: {output_path}")

if __name__ == "__main__":
    input_file = "data/reaction_dataset/C_N_Coupling.jsonl"
    output_file = "data/reaction_dataset/C_N_Coupling_Recommendation_V2.jsonl"
    
    # Process all reactions
    process_recommendation_dataset(input_file, output_file, limit=-1)
