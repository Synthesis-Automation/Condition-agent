#!/usr/bin/env python3
"""
Convert HTE CSV to a canonical, simplified format using the project's taxonomy.
Renames AREA_TOTAL_REDUCED to yield and maps reactant types to high-precision motifs.
"""

import csv
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# Mappings for Reaction Types
_REACTION_MAPPING = {
    "ch-activation": "c_h_activation",
    "suzuki": "suzuki_miyaura",
    "c_n_coupling": "cn_coupling",
    "amide_formation": "amide_formation",
    "snar": "snar",
    "heck": "heck",
    "sonogashira": "sonogashira",
    "negishi": "negishi",
    "stille": "stille",
    "cs-coupling": "cs_coupling",
    "co-coupling": "co_coupling",
    "condensation": "condensation",
    "cyclization": "cyclization",
    "borylation-miyaura": "miyaura_borylation",
    "borylation-c-h": "c_h_borylation",
    "reduction": "reduction",
    "oxidation": "oxidation",
    "alkylation": "alkylation",
    "halogenation-chlorination": "chlorination",
    "cyanation": "cyanation",
    "wittig": "wittig",
}

# Mappings for Reactant Types (High-Precision)
_COMPOUND_MAPPING = {
    # Halides
    "arbr": "Ar-Br",
    "arcl": "Ar-Cl",
    "ari": "Ar-I",
    "arf": "Ar-F",
    "aroso2r": "Ar-OSO2R",
    "arots": "Ar-OTs",
    "arotf": "Ar-OTf",
    "aroms": "Ar-OMs",
    "alkyl-oso2r": "R-OSO2R",
    "alkyl-cl": "R-Cl",
    "alkyl-br": "R-Br",
    "alkyl-i": "R-I",
    "alkyl-f": "R-F",
    "alkyl-x": "R-X",
    "alkene-i": "Alkenyl-I",
    "alkene-br": "Alkenyl-Br",
    
    # Boron
    "arb(or)2": "Ar-B(OR)2",
    "arb(oh)2": "Ar-B(OH)2",
    "arbf3k": "Ar-BF3K",
    "alkeneb(or)2": "Alkenyl-B(OR)2",
    "alkyl-b(oh)2": "RCH2-B(OH)2",
    "alkyl-b(or)2": "RCH2-B(OR)2",
    
    # Nitrogen
    "arnh2": "Ar-NH2",
    "arnhr": "Ar-NHR",
    "ar2nh": "Ar-NR2",
    "arom-nh": "AromN-H",
    "arom.nh": "AromN-H",
    "aromnh": "AromN-H",
    "ar-nh": "AromN-H",
    "rnh2": "R-NH2",
    "r2nh": "R-NHR",
    "r3n": "R-NR2",
    "rnh2-a-branch": "R2CH-NH2",
    "r2nh-a-branch": "R2CH-NHR",
    "rnh2a-branch": "R2CH-NH2",
    "r2nha-branch": "R2CH-NHR",
    "r2cnr": "R2CH-NHR",
    "rnhnh2": "R-Hydrazine",
    "rconh2": "R-CONH2",
    "rconhr": "R-CONHR",
    "lactam": "Any-CONHR",
    "carbamate": "R-Carbamate",
    "roconr2": "R-Carbamate",
    "roc(o)nr2": "R-Carbamate",
    "urea": "Any-Urea",
    "thiourea": "R-Thiourea",
    
    # Oxygen/Sulfur
    "aroh": "Ar-OH",
    "roh-primary": "RCH2-OH",
    "roh-a-branch": "R2CH-OH",
    "alkyl-ohprimary": "RCH2-OH",
    "alkyl-oha-branch": "R2CH-OH",
    "alkyl-oh": "R-OH",
    "ror": "R-OR",
    "rco2h": "R-CO2H",
    "rco2horm": "R-CO2H",
    "rco2r": "R-Ester",
    "rcho": "R-CHO",
    "r2co": "R-COR",
    "rsh": "R-SH",
    "rso2cl": "R-SO2Cl",
    "rsnr3": "R-Sn",
    "rocssr": "R-OCSSR",
    
    # Hydrocarbons
    "arh": "Ar-H",
    "ar-h": "Ar-H",
    "alkyl-h": "R-H",
    "alkyne": "R-CCH",
    "alkene": "Any-Alkene",
    "enol-ether": "Any-OR",
    
    # Others
    "alkyl-m": "R-Mg",
    "rcn": "R-CN",
    "rch2pph3x": "RCH2-PR3+",
}

def normalize_token(token: str) -> str:
    if not token:
        return ""
    token = token.strip().lower()
    # Remove common noise (stars, spaces, dots)
    token = re.sub(r"[\*\s\.]+", "", token)
    return _COMPOUND_MAPPING.get(token, token)

def split_and_map(value: str) -> List[str]:
    if not value:
        return []
    # Split by comma or semicolon
    tokens = re.split(r"[,;]+", value)
    results = []
    for t in tokens:
        mapped = normalize_token(t)
        if mapped:
            results.append(mapped)
    return sorted(list(set(results)))

def convert_hte(input_path: str, output_path: str):
    input_file = Path(input_path)
    output_file = Path(output_path)
    
    if not input_file.exists():
        print(f"Error: {input_path} not found.")
        return

    print(f"Converting {input_path} to {output_path}...")
    
    rows_processed = 0
    with input_file.open("r", encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in)
        
        # Define new simplified headers
        headers = [
            "reaction_type", "yield", "z_score", 
            "reactant_1", "orig_reactant_1",
            "reactant_2", "orig_reactant_2",
            "catalyst", "ligand", "base", "solvent", "additive"
        ]
        
        with output_file.open("w", encoding="utf-8", newline="") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=headers)
            writer.writeheader()
            
            for row in reader:
                # Skip rows containing "protectinggroup" (case-insensitive, ignoring spaces)
                row_str = "".join(str(v) for v in row.values()).lower().replace(" ", "")
                if "protectinggroup" in row_str:
                    continue

                # 1. Map Reaction Type
                raw_rxn = (row.get("Reaction_Type_Standardized") or "").strip().lower()
                rxn_type = _REACTION_MAPPING.get(raw_rxn, raw_rxn)
                
                # 2. Map Reactants
                # Try to get reactant 1 from A label, then A type
                r_a_tokens = split_and_map(row.get("Reactant_A"))
                orig_r1 = row.get("Reactant_A")
                if not r_a_tokens:
                    r_a_tokens = split_and_map(row.get("Reactant_A_Type"))
                    orig_r1 = row.get("Reactant_A_Type")
                
                # Try to get reactant 2 from B label, then B type
                r_b_tokens = split_and_map(row.get("Reactant_B"))
                orig_r2 = row.get("Reactant_B")
                if not r_b_tokens:
                    r_b_tokens = split_and_map(row.get("Reactant_B_Type"))
                    orig_r2 = row.get("Reactant_B_Type")
                
                # Assign to slots
                reactant_1 = r_a_tokens[0] if r_a_tokens else ""
                reactant_2 = r_b_tokens[0] if r_b_tokens else ""
                
                # Special case: if they are the same, but one side has an alternative in the Type column
                if reactant_1 == reactant_2 and reactant_1:
                    alt_a = split_and_map(row.get("Reactant_A_Type"))
                    alt_b = split_and_map(row.get("Reactant_B_Type"))
                    if alt_a and alt_a[0] != reactant_1:
                        reactant_1 = alt_a[0]
                        orig_r1 = row.get("Reactant_A_Type")
                    elif alt_b and alt_b[0] != reactant_2:
                        reactant_2 = alt_b[0]
                        orig_r2 = row.get("Reactant_B_Type")
                
                # 3. Clean Conditions
                catalyst = (row.get("Catalyst") or "").strip()
                ligand = (row.get("Ligand") or "").strip()
                base = (row.get("Base") or "").strip()
                
                # Combine solvents
                solvents = [s for s in [row.get("Solvent"), row.get("Secondary Solvent")] if s and s.strip()]
                solvent = ", ".join(solvents)
                
                # Combine additives and coupling reagents
                additives = [a for a in [row.get("Additive"), row.get("Coupling Reagent")] if a and a.strip()]
                additive = ", ".join(additives)
                
                # 4. Metrics
                yield_val = row.get("AREA_TOTAL_REDUCED")
                z_score = row.get("z-Score")
                
                writer.writerow({
                    "reaction_type": rxn_type,
                    "yield": yield_val,
                    "z_score": z_score,
                    "reactant_1": reactant_1,
                    "orig_reactant_1": orig_r1,
                    "reactant_2": reactant_2,
                    "orig_reactant_2": orig_r2,
                    "catalyst": catalyst,
                    "ligand": ligand,
                    "base": base,
                    "solvent": solvent,
                    "additive": additive
                })
                
                rows_processed += 1
                if rows_processed % 10000 == 0:
                    print(f"Processed {rows_processed} rows...")

    print(f"Done! Converted {rows_processed} rows.")

if __name__ == "__main__":
    convert_hte("data/HTE_db/HTE_0.csv", "data/HTE_db/HTE_canonical.csv")
