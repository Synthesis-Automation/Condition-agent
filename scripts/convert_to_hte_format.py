import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Tuple, Iterable, Optional
import numpy as np
from functools import lru_cache

# Import chemtools components
from chemtools.featurizers.structural import featurize_molecule

@lru_cache(maxsize=10000)
def cached_featurize(smiles: str):
    """Cache featurization results to speed up processing of repeated reactants."""
    return featurize_molecule(smiles)

def _dedupe_list(values: Iterable[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for value in values:
        if not value:
            continue
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out

def _reactant_key(values: Iterable[Optional[str]]) -> str:
    """Generate a standardized key for reactant combinations."""
    items = _dedupe_list([str(v).strip() for v in values if v])
    return "|".join(sorted(items))

def extract_reagents(record: Dict[str, Any]) -> Dict[str, str]:
    """Extract and categorize reagents from the reaction record."""
    reagents = record.get("reagents", [])
    solvents = record.get("solvents", [])
    
    bases = [r["name"] for r in reagents if r.get("role") == "BASE"]
    additives = [r["name"] for r in reagents if r.get("role") == "ADDITIVE"]
    coupling_reagents = [r["name"] for r in reagents if r.get("role") == "COUPLING_REAGENT"]
    
    solvent_names = [s["name"] for s in solvents]
    
    # Catalyst and Ligand extraction heuristics
    catalytic_system = record.get("catalytic_system", "")
    condition_core = record.get("condition_core", "")
    
    catalyst = ""
    ligand = ""
    
    if catalytic_system:
        parts = [p.strip() for p in catalytic_system.split(",")]
        if len(parts) >= 1:
            catalyst = parts[0]
        if len(parts) >= 2:
            ligand = parts[1]
    elif condition_core and "/" in condition_core:
        parts = [p.strip() for p in condition_core.split("/")]
        # Heuristic: if first part looks like a metal
        metals = ["Pd", "Ni", "Cu", "Ir", "Rh", "Ru", "Pt", "Au", "Ag", "Fe", "Co", "Zn"]
        if any(m in parts[0] for m in metals):
            catalyst = parts[0]
            if len(parts) >= 2:
                ligand = parts[1]
    
    return {
        "Catalyst": catalyst,
        "Ligand": ligand,
        "Base": " / ".join(bases),
        "Solvent": " / ".join(solvent_names),
        "Additive": " / ".join(additives),
        "Coupling Reagent": " / ".join(coupling_reagents),
        "Secondary Solvent": ""
    }

def process_reaction_dataset(input_path: str, output_path: str):
    """Convert reaction dataset to HTE-canonical CSV format."""
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    rows = []
    print(f"Reading {input_path}...")
    if not input_path.exists():
        print(f"Error: Input file {input_path} not found.")
        return

    with open(input_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    total = len(lines)
    print(f"Processing {total} reactions...")
    
    for i, line in enumerate(lines):
        if i % 500 == 0:
            print(f"Progress: {i}/{total} ({(i/total)*100:.1f}%)")
            
        try:
            record = json.loads(line)
        except:
            continue
            
        smiles = record.get("smiles", "")
        if ">>" not in smiles:
            continue
            
        reactants_part = smiles.split(">>")[0]
        reactants = reactants_part.split(".")
        
        motif_ids = []
        reactant_data = []
        
        for r_smiles in reactants:
            try:
                # Use cached featurization for speed
                analysis = cached_featurize(r_smiles)
                motifs = analysis.get("motifs", [])
                
                r_motifs = _dedupe_list([m.get("compound_id", "") for m in motifs if m.get("compound_id")])
                r_cat = motifs[0].get("category", "Unknown") if motifs else "Unknown"
                
                reactant_data.append({
                    "motifs": r_motifs,
                    "category": r_cat
                })
                motif_ids.extend(r_motifs)
            except Exception as e:
                # Skip invalid SMILES or featurization errors
                continue
            
        if not reactant_data:
            continue
            
        # Map to A and B slots (standard for HTE recommender)
        type_a = ",".join(reactant_data[0]["motifs"]) if len(reactant_data) > 0 else ""
        
        type_b = ",".join(reactant_data[1]["motifs"]) if len(reactant_data) > 1 else ""
        
        reagents = extract_reagents(record)
        
        row = {
            "Reaction_Type_Standardized": record.get("type", "Unknown"),
            "Reactant_A_Type": type_a,
            "Reactant_B_Type": type_b,
            "Reactant_Types_Key": _reactant_key(motif_ids),
            "yield": record.get("yield", 0.0),
            "smiles": smiles,
            **reagents
        }
        rows.append(row)
        
    df = pd.DataFrame(rows)
    
    # Ensure yield is numeric
    df["yield"] = pd.to_numeric(df["yield"], errors='coerce').fillna(0.0)
    
    # Calculate Pseudo Z-Score based on motif groups
    print("Calculating Pseudo Z-Scores...")
    df["z-Score"] = 0.0
    df["z-Score"] = df["z-Score"].astype(float)
    
    for key, group in df.groupby("Reactant_Types_Key"):
        if len(group) > 1:
            yields = group["yield"]
            mean = yields.mean()
            std = yields.std()
            if std > 0:
                df.loc[group.index, "z-Score"] = (yields - mean) / std
            else:
                # If all yields are same, Z-score is 0
                df.loc[group.index, "z-Score"] = 0.0
        else:
            # Single experiment in group
            df.loc[group.index, "z-Score"] = 0.0
            
    df.to_csv(output_path, index=False)
    print(f"Successfully saved {len(df)} reactions to {output_path}")

if __name__ == "__main__":
    # Process the C-N Coupling dataset
    process_reaction_dataset(
        "data/reaction_dataset/C_N_Coupling.jsonl",
        "data/HTE_db/C_N_Coupling_canonical.csv"
    )
