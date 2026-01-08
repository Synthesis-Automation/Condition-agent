import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Tuple, Iterable, Optional
import numpy as np
from functools import lru_cache
import argparse
import sys
import os

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

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

def process_reaction_dataset(input_path: str, output_path: str, drop_no_catalyst: bool = True):
    """Convert reaction dataset to HTE-canonical CSV format with extra metadata."""
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
                current_r_motifs = []
                for m in motifs:
                    cid = m.get("compound_id", "")
                    if cid:
                        current_r_motifs.append(cid)
                        # Include alternatives to ensure information is preserved across perspectives
                        for alt_id in m.get("alt_compound_ids", []):
                            current_r_motifs.append(alt_id)
                
                reactant_data.append({
                    "motifs": _dedupe_list(current_r_motifs)
                })
                motif_ids.extend(current_r_motifs)
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
            for m in p_motifs:
                cid = m.get("compound_id", "")
                if cid:
                    product_motifs.append(cid)
                    for alt_id in m.get("alt_compound_ids", []):
                        product_motifs.append(alt_id)
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
        
        # Picking ONE motif for each reactant and product for the reaction_key
        # For each reactant, pick the highest priority motif that actually reacted
        primary_reacted_motifs = []
        for r_info in reactant_data:
            r_motifs = r_info["motifs"]
            # Find which of these are in the reacted_set
            reacted_here = [m for m in r_motifs if m in reacted_set]
            if reacted_here:
                # Pick the first one (highest priority)
                primary_reacted_motifs.append(reacted_here[0])
        
        primary_reacted_str = _reactant_key(primary_reacted_motifs) or "None"
        
        # For the product, pick the highest priority motif that was formed
        formed_here = [m for m in _dedupe_list(product_motifs) if m in formed_set]
        primary_formed_str = formed_here[0] if formed_here else "None"
        
        # Spectators (unchanged motifs)
        spectators_str = _reactant_key(list(spectators_set)) or "None"

        # Primary motif per reactant (fallback to first motif if none reacted)
        primary_reactant_motifs = []
        for r_info in reactant_data:
            r_motifs = r_info["motifs"]
            reacted_here = [m for m in r_motifs if m in reacted_set]
            if reacted_here:
                primary_reactant_motifs.append(reacted_here[0])
            elif r_motifs:
                primary_reactant_motifs.append(r_motifs[0])
            else:
                primary_reactant_motifs.append("")
        
        # Combined reaction key: A|B -> P || Unchanged
        # This replaces Reaction_Type_Standardized and Reactant_Types_Key
        reaction_key = f"{primary_reacted_str} -> {primary_formed_str} || {spectators_str}"
            
        # All motifs for reference
        all_reactant_motifs_str = "|".join(sorted(_dedupe_list(motif_ids)))
        all_product_motifs_str = "|".join(sorted(_dedupe_list(product_motifs)))

        reacted_motifs_str = _reactant_key(list(reacted_set)) or "None"
        formed_motifs_str = _reactant_key(list(formed_set)) or "None"
        
        # Map to A and B slots (canonical for HTE recommender)
        type_a = primary_reactant_motifs[0] if len(primary_reactant_motifs) > 0 else ""
        type_b = primary_reactant_motifs[1] if len(primary_reactant_motifs) > 1 else ""
        type_p = ",".join(_dedupe_list(product_motifs))
        
        reagents = extract_reagents(record)
        if drop_no_catalyst and not reagents.get("Catalyst"):
            continue
            
        row = {
            # HTE canonical columns (lowercase with underscores)
            "reaction_type": reaction_key,
            "yield": record.get("yield", 0.0),
            "z_score": 0.0,
            "reactant_1": type_a,
            "reactant_2": type_b,
            "catalyst": reagents.get("Catalyst", ""),
            "ligand": reagents.get("Ligand", ""),
            "base": reagents.get("Base", ""),
            "solvent": reagents.get("Solvent", ""),
            "additive": reagents.get("Additive", ""),

            # Extra metadata for downstream analysis
            "reaction_smiles": smiles,
            "smiles": smiles,
            "reactant_smiles": reactants_part,
            "product_smiles": products_part,
            "reaction_key": reaction_key,
            "reactant_types_key": _reactant_key([type_a, type_b]),
            "product_type": type_p,
            "reactant_a_all_motifs": ",".join(reactant_data[0]["motifs"]) if len(reactant_data) > 0 else "",
            "reactant_b_all_motifs": ",".join(reactant_data[1]["motifs"]) if len(reactant_data) > 1 else "",
            "all_reactant_motifs": all_reactant_motifs_str,
            "all_product_motifs": all_product_motifs_str,
            "reacted_motifs": reacted_motifs_str,
            "formed_motifs": formed_motifs_str,
            "spectator_motifs": spectators_str,
            "is_intramolecular": len(reactants) == 1,
            "coupling_reagent": reagents.get("Coupling Reagent", ""),
            "secondary_solvent": reagents.get("Secondary Solvent", ""),
            "catalytic_system": record.get("catalytic_system", ""),
            "condition_core": record.get("condition_core", ""),
            "raw_reagents": json.dumps(record.get("reagents", []), ensure_ascii=True),
            "raw_solvents": json.dumps(record.get("solvents", []), ensure_ascii=True),
        }
        rows.append(row)
        
    df = pd.DataFrame(rows)
    
    # Ensure yield is numeric
    df["yield"] = pd.to_numeric(df["yield"], errors='coerce').fillna(0.0)
    
    # Calculate Pseudo Z-Score based on motif groups
    print("Calculating Pseudo Z-Scores...")
    df["z_score"] = 0.0
    df["z_score"] = df["z_score"].astype(float)
    
    for key, group in df.groupby("reaction_type"):
        if len(group) > 1:
            yields = group["yield"]
            mean = yields.mean()
            std = yields.std()
            if std > 0:
                df.loc[group.index, "z_score"] = (yields - mean) / std
            else:
                # If all yields are same, Z-score is 0
                df.loc[group.index, "z_score"] = 0.0
        else:
            # Single experiment in group
            df.loc[group.index, "z_score"] = 0.0

    canonical_cols = [
        "reaction_type", "yield", "z_score", "reactant_1", "reactant_2",
        "catalyst", "ligand", "base", "solvent", "additive"
    ]
    extra_cols = [col for col in df.columns if col not in canonical_cols]
    df = df[canonical_cols + extra_cols]

    df.to_csv(output_path, index=False)
    print(f"Successfully saved {len(df)} reactions to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert reaction dataset to HTE-canonical CSV format.")
    parser.add_argument("--dataset", "-d", help="Name of the dataset (e.g., C_N_Coupling). If provided, processes 'data/reaction_dataset/{dataset}.jsonl' to 'data/HTE_db/{dataset}_canonical.csv'.")
    parser.add_argument("--input", "-i", help="Direct path to input JSONL file.")
    parser.add_argument("--output", "-o", help="Direct path to output CSV file.")
    
    args = parser.parse_args()
    
    if args.dataset:
        input_file = f"data/reaction_dataset/{args.dataset}.jsonl"
        output_file = f"data/HTE_db/{args.dataset}_canonical.csv"
        process_reaction_dataset(input_file, output_file)
    elif args.input and args.output:
        process_reaction_dataset(args.input, args.output)
    else:
        # Default: Process known core datasets
        datasets = ["C_N_Coupling", "C_O_Coupling", "C_S_Coupling"]
        for ds in datasets:
            input_file = f"data/reaction_dataset/{ds}.jsonl"
            output_file = f"data/HTE_db/{ds}_canonical.csv"
            process_reaction_dataset(input_file, output_file)
