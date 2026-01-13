import csv
import json
import re
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Tuple, Iterable, Optional
import numpy as np
from functools import lru_cache
import argparse
import sys
import os
from collections import Counter, defaultdict

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Import chemtools components
from chemtools.featurizers.structural import featurize_molecule
from chemtools.smiles import normalize

@lru_cache(maxsize=10000)
def cached_featurize(smiles: str):
    """Cache featurization results to speed up processing of repeated reactants."""
    # Use substituent-level filtering to avoid redundant labels (e.g., Ar-NR2 and RCH2-NR2)
    # also ensures that Aromatic scaffolds win over Aliphatic ones due to updated priorities.
    options = {"motif_site_filter": "substituent"}
    return featurize_molecule(smiles, options=options)

CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")
CAS_INLINE_PATTERN = re.compile(r"\b\d{2,7}-\d{2}-\d\b")

def _normalize_smiles(smiles: str) -> Optional[str]:
    if not smiles:
        return None
    info = normalize(smiles)
    return info.get("smiles_norm") or info.get("largest_smiles") or None

def _split_reagent_names(value: str) -> List[str]:
    if not value:
        return []
    parts = []
    for token in value.replace(";", "/").split("/"):
        for item in token.split(","):
            name = item.strip()
            if name:
                parts.append(name)
    return parts

def _normalize_reagent_text(value: str) -> str:
    return " ".join(str(value).strip().lower().split())

@lru_cache(maxsize=4)
def _load_reagent_csv(path_str: str) -> Dict[str, Any]:
    path = Path(path_str)
    data = {
        "name_to_smiles": defaultdict(set),
        "cas_to_smiles": defaultdict(set),
        "cas_set": set(),
    }
    if not path.exists():
        return data
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = _normalize_reagent_text(row.get("name", ""))
            abbr = _normalize_reagent_text(row.get("abbreviation", ""))
            cas = (row.get("cas") or "").strip()
            smile = (row.get("smiles") or row.get("smile") or "").strip()
            if cas:
                data["cas_set"].add(cas)
            if smile:
                if name:
                    data["name_to_smiles"][name].add(smile)
                if abbr:
                    data["name_to_smiles"][abbr].add(smile)
                if cas:
                    data["cas_to_smiles"][cas].add(smile)
    return data

def _collect_reagent_smiles(
    record: Dict[str, Any],
    reagent_db: Dict[str, Any],
    unknown_cas: Optional[set[str]] = None,
) -> set[str]:
    if not reagent_db:
        return set()

    def normalize_entries(values: Any) -> List[Dict[str, Any]]:
        if values is None:
            return []
        if isinstance(values, dict):
            values = [values]
        entries = []
        for entry in values:
            if isinstance(entry, dict):
                entries.append(entry)
            elif isinstance(entry, str):
                entries.append({"name": entry})
        return entries

    reagent_names: List[str] = []
    reagent_cas: List[str] = []
    for entry in normalize_entries(record.get("reagents")):
        name = entry.get("name")
        if name:
            reagent_names.extend(_split_reagent_names(str(name)))
        cas = entry.get("cas")
        if cas:
            reagent_cas.append(str(cas).strip())
    for entry in normalize_entries(record.get("solvents")):
        name = entry.get("name")
        if name:
            reagent_names.extend(_split_reagent_names(str(name)))
        cas = entry.get("cas")
        if cas:
            reagent_cas.append(str(cas).strip())

    catalytic_system = record.get("catalytic_system") or ""
    condition_core = record.get("condition_core") or ""
    if catalytic_system:
        reagent_names.extend(_split_reagent_names(str(catalytic_system)))
    if condition_core:
        reagent_names.extend(_split_reagent_names(str(condition_core)))

    reagent_smiles = set()
    cas_set = reagent_db.get("cas_set", set())
    cas_to_smiles = reagent_db.get("cas_to_smiles", {})
    name_to_smiles = reagent_db.get("name_to_smiles", {})

    unknown = unknown_cas if unknown_cas is not None else set()

    for cas in _dedupe_list(reagent_cas):
        if not cas or cas.startswith("$") or not CAS_PATTERN.match(cas):
            continue
        if cas not in cas_set:
            unknown.add(cas)
        for smiles in cas_to_smiles.get(cas, []):
            normalized = _normalize_smiles(smiles) if smiles else None
            if normalized:
                reagent_smiles.add(normalized)

    for name in _dedupe_list(reagent_names):
        key = _normalize_reagent_text(name)
        for match in CAS_INLINE_PATTERN.findall(key):
            reagent_cas.append(match)
        for smiles in name_to_smiles.get(key, []):
            normalized = _normalize_smiles(smiles) if smiles else None
            if normalized:
                reagent_smiles.add(normalized)
    return reagent_smiles

def _write_new_reagents(path: Path, cas_values: set[str]) -> None:
    if not cas_values:
        return
    existing: set[str] = set()
    if path.exists():
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                cas = (row.get("cas") or "").strip()
                if cas:
                    existing.add(cas)
    new_only = {cas.strip() for cas in cas_values if cas.strip()} - existing
    if not new_only:
        return
    merged = sorted(existing | new_only)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cas"])
        writer.writeheader()
        for cas in merged:
            writer.writerow({"cas": cas})

@lru_cache(maxsize=1)
def _load_scaffold_motif_ids() -> set[str]:
    path = PROJECT_ROOT / "chemtools" / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    motifs = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = str(entry.get("id") or "").strip()
        if motif_id:
            motifs.add(motif_id)
    return motifs

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

def _motif_group_ids(values: Iterable[str]) -> set[str]:
    group_ids: set[str] = set()
    for value in values:
        token = str(value).strip()
        if not token:
            continue
        if "-" in token:
            group_id = token.split("-")[-1]
        else:
            group_id = token
        group_id = group_id.strip()
        if group_id:
            group_ids.add(group_id)
    return group_ids

_SPECTATOR_GROUP_STOPLIST = {
    "Ar",
    "R",
    "Any",
    "Alkyl",
    "Alkenyl",
    "Alkynyl",
    "H",
}

def _group_id_from_motif_id(motif_id: str) -> str:
    text = str(motif_id).strip()
    if not text:
        return ""
    if "-" in text:
        return text.split("-")[-1].strip()
    return text

def _collect_spectator_groups(
    reactant_data: List[Dict[str, Any]],
    spectators_set: set[str],
) -> List[str]:
    scaffold_ids = _load_scaffold_motif_ids()
    seen = set()
    groups: List[str] = []

    for r_info in reactant_data:
        motifs = r_info.get("motifs") or []
        for motif_id in motifs:
            if motif_id not in spectators_set:
                continue
            group_id = _group_id_from_motif_id(motif_id)
            if not group_id or group_id in _SPECTATOR_GROUP_STOPLIST:
                continue
            if group_id in seen:
                continue
            seen.add(group_id)
            groups.append(group_id)

    for motif_id in spectators_set:
        if motif_id in scaffold_ids and motif_id not in seen:
            seen.add(motif_id)
            groups.append(motif_id)

    return groups

def extract_reagents(record: Dict[str, Any]) -> Dict[str, str]:
    """Extract and categorize reagents from the reaction record."""
    def normalize_entries(values: Any) -> List[Dict[str, Any]]:
        if values is None:
            return []
        if isinstance(values, dict):
            values = [values]
        entries = []
        for entry in values:
            if isinstance(entry, dict):
                entries.append(entry)
            elif isinstance(entry, str):
                entries.append({"name": entry})
        return entries

    reagents = normalize_entries(record.get("reagents"))
    solvents = normalize_entries(record.get("solvents"))
    
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

def process_reaction_dataset(
    input_path: str,
    output_path: str,
    drop_no_catalyst: bool = True,
    drop_reagent_reactants: bool = True,
    reagent_csv_path: Optional[str | Path] = None,
    new_reagents_path: Optional[str | Path] = None,
):
    """Convert reaction dataset to a minimal HTE recommender CSV."""
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_label = input_path.stem

    reagent_csv = Path(reagent_csv_path) if reagent_csv_path else (
        PROJECT_ROOT / "data" / "reagent_db" / "reagents.csv"
    )
    reagent_db = _load_reagent_csv(str(reagent_csv))
    unknown_cas: set[str] = set()
    new_reagents_csv = Path(new_reagents_path) if new_reagents_path else (
        PROJECT_ROOT / "data" / "reagent_db" / "new_reagents.csv"
    )
    
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
        reagent_smiles = _collect_reagent_smiles(record, reagent_db, unknown_cas)
        if drop_reagent_reactants and reagent_smiles:
            filtered = []
            for r_smiles in reactants:
                norm = _normalize_smiles(r_smiles)
                if norm and norm in reagent_smiles:
                    continue
                filtered.append(r_smiles)
            if filtered:
                reactants = filtered
        
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
                    "motifs": _dedupe_list(current_r_motifs),
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
            
        # Map to A and B slots (canonical for HTE recommender)
        type_a = primary_reactant_motifs[0] if len(primary_reactant_motifs) > 0 else ""
        type_b = primary_reactant_motifs[1] if len(primary_reactant_motifs) > 1 else ""
        
        reagents = extract_reagents(record)
        if drop_no_catalyst and not reagents.get("Catalyst"):
            continue

        spectator_groups = _collect_spectator_groups(reactant_data, spectators_set)
            
        row = {
            # HTE canonical columns (lowercase with underscores)
            "reaction_id": source_label,
            "yield": record.get("yield", 0.0),
            "z_score": 0.0,
            "reactant_1": type_a,
            "reactant_2": type_b,
            "catalyst": reagents.get("Catalyst", ""),
            "ligand": reagents.get("Ligand", ""),
            "base": reagents.get("Base", ""),
            "solvent": reagents.get("Solvent", ""),
            "additive": reagents.get("Additive", ""),
            "Secondary Solvent": reagents.get("Secondary Solvent", ""),
            "Coupling Reagent": reagents.get("Coupling Reagent", ""),
            "Is_Intramolecular": len(reactants) == 1,
            "reaction_smiles": smiles,
            "spectator_groups": " / ".join(spectator_groups),
            "_reaction_key": reaction_key,
        }
        rows.append(row)
        
    df = pd.DataFrame(rows)
    
    # Ensure yield is numeric
    df["yield"] = pd.to_numeric(df["yield"], errors='coerce').fillna(0.0)
    
    # Calculate Pseudo Z-Score based on motif groups
    print("Calculating Pseudo Z-Scores...")
    df["z_score"] = 0.0
    df["z_score"] = df["z_score"].astype(float)
    
    if "_reaction_key" in df.columns:
        group_key = "_reaction_key"
    else:
        group_key = "reaction_id"

    for key, group in df.groupby(group_key):
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

    df["z_score"] = df["z_score"].round(2)

    canonical_cols = [
        "reaction_id", "yield", "z_score", "reactant_1", "reactant_2",
        "catalyst", "ligand", "base", "solvent", "additive",
        "Secondary Solvent", "Coupling Reagent", "Is_Intramolecular",
        "reaction_smiles", "spectator_groups",
    ]
    df = df[canonical_cols]

    df.to_csv(output_path, index=False)
    print(f"Successfully saved {len(df)} reactions to {output_path}")
    if unknown_cas:
        _write_new_reagents(new_reagents_csv, unknown_cas)
        print(f"Saved {len(unknown_cas)} new reagent CAS entries to {new_reagents_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert reaction dataset to HTE-canonical CSV format.")
    parser.add_argument("--dataset", "-d", help="Name of the dataset (e.g., C_N_Coupling). If provided, processes 'data/reaction_dataset/{dataset}.jsonl' to 'data/HTE_db/{dataset}_canonical.csv'.")
    parser.add_argument("--input", "-i", help="Direct path to input JSONL file.")
    parser.add_argument("--output", "-o", help="Direct path to output CSV file.")
    parser.add_argument(
        "--keep-reagent-reactants",
        action="store_true",
        help="Keep reagent/solvent molecules in the reactant columns (default: drop).",
    )
    parser.add_argument(
        "--reagent-csv",
        help="Path to reagent registry CSV (default: data/reagent_db/reagents.csv).",
    )
    parser.add_argument(
        "--new-reagents",
        help="Path to write missing reagent CAS list (default: data/reagent_db/new_reagents.csv).",
    )
    
    args = parser.parse_args()
    
    drop_reagent_reactants = not args.keep_reagent_reactants

    if args.dataset:
        input_file = f"data/reaction_dataset/{args.dataset}.jsonl"
        output_file = f"data/HTE_db/{args.dataset}_canonical.csv"
        process_reaction_dataset(
            input_file,
            output_file,
            drop_reagent_reactants=drop_reagent_reactants,
            reagent_csv_path=args.reagent_csv,
            new_reagents_path=args.new_reagents,
        )
    elif args.input and args.output:
        process_reaction_dataset(
            args.input,
            args.output,
            drop_reagent_reactants=drop_reagent_reactants,
            reagent_csv_path=args.reagent_csv,
            new_reagents_path=args.new_reagents,
        )
    else:
        # Default: Process known core datasets
        datasets = ["C_N_Coupling", "C_O_Coupling", "C_S_Coupling"]
        for ds in datasets:
            input_file = f"data/reaction_dataset/{ds}.jsonl"
            output_file = f"data/HTE_db/{ds}_canonical.csv"
            process_reaction_dataset(
                input_file,
                output_file,
                drop_reagent_reactants=drop_reagent_reactants,
                reagent_csv_path=args.reagent_csv,
                new_reagents_path=args.new_reagents,
            )
