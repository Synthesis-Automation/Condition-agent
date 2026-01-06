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
    # Use substituent-level filtering to avoid redundant labels (e.g., Ar-NR2 and RCH2-NR2)
    # also ensures that Aromatic scaffolds win over Aliphatic ones due to updated priorities.
    options = {"motif_site_filter": "substituent"}
    return featurize_molecule(smiles, options=options)

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
            
        reactants_part, products_part = smiles.split(">>")
        reactants = reactants_part.split(".")
        
        motif_ids = []
        reactant_data = []
        
        for r_smiles in reactants:
            try:
                # Use cached featurization for speed
                analysis = cached_featurize(r_smiles)
                motifs = analysis.get("motifs", [])
                
                # Keep raw motifs for stoichiometry counts
                r_motifs = [m.get("compound_id", "") for m in motifs if m.get("compound_id")]
                
                reactant_data.append({
                    "motifs": _dedupe_list(r_motifs) # Dedupe for the per-reactant display
                })
                motif_ids.extend(r_motifs)
            except Exception as e:
                # Skip invalid SMILES or featurization errors
                continue
            
        if not reactant_data:
            continue

        # Product motif analysis
        product_motifs = []
        try:
            p_analysis = cached_featurize(products_part)
            p_motifs = p_analysis.get("motifs", [])
            product_motifs = [m.get("compound_id", "") for m in p_motifs if m.get("compound_id")]
        except:
            pass

        # Transformation analysis using counts to handle motifs that exist in both
        from collections import Counter
        r_counts = Counter(motif_ids)
        p_counts = Counter(product_motifs)
        
        # Heuristic: certain common heterocycles often appear as reagents/bases/solvents
        # and may be missing from the product SMILES even if they don't react.
        # We treat them as persistent spectators if they disappear.
        PERSISTENT_HETEROCYCLES = {"Pyridine", "Quinoline", "Isoquinoline", "Pyrimidine", "Pyrazine", "Pyridazine"}
        for h_id in PERSISTENT_HETEROCYCLES:
            if r_counts[h_id] > 0 and p_counts[h_id] < r_counts[h_id]:
                p_counts[h_id] = r_counts[h_id]

        reacted_set = set()
        formed_set = set()
        spectators_set = set()
        
        all_motifs = set(r_counts.keys()) | set(p_counts.keys())
        for m in all_motifs:
            rc = r_counts.get(m, 0)
            pc = p_counts.get(m, 0)
            
            if pc > rc:
                formed_set.add(m)
                if rc > 0:
                    spectators_set.add(m)
            elif pc < rc:
                reacted_set.add(m)
                if pc > 0:
                    spectators_set.add(m)
            else:
                if rc > 0:
                    spectators_set.add(m)
        
        # Construct the new informative key
        # Format: [Reacted] -> [Formed] || [Spectators]
        reacted_str = _reactant_key(list(reacted_set)) or "None"
        formed_str = _reactant_key(list(formed_set)) or "None"
        spectators_str = _reactant_key(list(spectators_set)) or "None"
        
        transformation_key = f"{reacted_str} -> {formed_str} || {spectators_str}"
            
        # Map to A and B slots (standard for HTE recommender)
        type_a = ",".join(reactant_data[0]["motifs"]) if len(reactant_data) > 0 else ""
        
        type_b = ",".join(reactant_data[1]["motifs"]) if len(reactant_data) > 1 else ""
        
        reagents = extract_reagents(record)
        
        row = {
            "Reaction_Type_Standardized": transformation_key,
            "Reactant_A_Type": type_a,
            "Reactant_B_Type": type_b,
            "Reactant_Types_Key": _reactant_key([type_a, type_b]),
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
