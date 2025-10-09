#!/usr/bin/env python3
"""
Fix Rule Database Reagents
===========================

This script ensures that all ligands, catalysts, and other reagents mentioned in
token_signature are properly included in the reagents array for each entry.
"""

import json
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def fix_suzuki_reagents(db_path: str):
    """Fix missing reagents in Suzuki database."""
    print(f"Processing: {db_path}")
    
    with open(db_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    ligand_mapping = {
        'SPhos': {'role': 'ligand', 'amount': '1.0–3.0%'},
        'XPhos': {'role': 'ligand', 'amount': '1.0–3.0%'},
        'PPh3': {'role': 'ligand', 'amount': '2.0–4.0%'},
        'tBuXPhos': {'role': 'ligand', 'amount': '1.0–3.0%'},
        'L95': {'role': 'ligand', 'amount': '1.0–3.0%'},
        'DPPF': {'role': 'ligand', 'amount': '2.0–5.0%'},
        'Amphos': {'role': 'ligand', 'amount': '1.0–3.0%'},
        'BrettPhos': {'role': 'ligand', 'amount': '1.0–3.0%'},
    }
    
    entries_fixed = 0
    
    for entry in data['entries']:
        token_sig = entry.get('token_signature', [])
        reagents = entry.get('reagents', [])
        reagent_names = {r.get('name', '') for r in reagents}
        
        # Check for missing ligands
        for lig_name, lig_info in ligand_mapping.items():
            if lig_name in token_sig and lig_name not in reagent_names:
                # Find position to insert (after metal source, before base)
                insert_pos = 1  # Default after metal source
                for idx, r in enumerate(reagents):
                    if r.get('role') in ['metal_source', 'metal_precursor', 'catalyst']:
                        insert_pos = idx + 1
                    elif r.get('role') == 'base':
                        break
                
                # Insert the ligand
                reagents.insert(insert_pos, {
                    'name': lig_name,
                    'role': lig_info['role'],
                    'amount': lig_info['amount']
                })
                print(f"  Fixed {entry['id']}: added {lig_name}")
                entries_fixed += 1
    
    # Save the fixed database
    with open(db_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"  Fixed {entries_fixed} entries")
    return entries_fixed


def fix_cn_coupling_reagents(db_path: str):
    """Fix missing reagents in C-N coupling databases."""
    print(f"Processing: {db_path}")
    
    with open(db_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    ligand_mapping = {
        # Cu ligands
        'L-proline': {'role': 'ligand', 'amount': '5–20%'},
        '1,10-phenanthroline': {'role': 'ligand', 'amount': '5–20%'},
        'DMEDA': {'role': 'ligand', 'amount': '10–20%'},
        'trans-1,2-cyclohexanediamine': {'role': 'ligand', 'amount': '10–20%'},
        '8-hydroxyquinoline': {'role': 'ligand', 'amount': '5–20%'},
        # Pd ligands
        'BINAP': {'role': 'ligand', 'amount': '2–5%'},
        'Xantphos': {'role': 'ligand', 'amount': '2–5%'},
        'XPhos': {'role': 'ligand', 'amount': '2–5%'},
        'tBuXPhos': {'role': 'ligand', 'amount': '2–5%'},
        'RuPhos': {'role': 'ligand', 'amount': '2–5%'},
        'SPhos': {'role': 'ligand', 'amount': '2–5%'},
        'BrettPhos': {'role': 'ligand', 'amount': '2–5%'},
        'dppf': {'role': 'ligand', 'amount': '2–5%'},
    }
    
    entries_fixed = 0
    
    for entry in data['entries']:
        token_sig = entry.get('token_signature', [])
        reagents = entry.get('reagents', [])
        reagent_names = {r.get('name', '') for r in reagents}
        
        # Check for missing ligands
        for lig_name, lig_info in ligand_mapping.items():
            if lig_name in token_sig and lig_name not in reagent_names:
                # Find position to insert
                insert_pos = 1
                for idx, r in enumerate(reagents):
                    if r.get('role') in ['metal_source', 'metal_precursor', 'catalyst']:
                        insert_pos = idx + 1
                    elif r.get('role') == 'base':
                        break
                
                # Insert the ligand
                reagents.insert(insert_pos, {
                    'name': lig_name,
                    'role': lig_info['role'],
                    'amount': lig_info['amount']
                })
                print(f"  Fixed {entry['id']}: added {lig_name}")
                entries_fixed += 1
    
    # Save the fixed database
    with open(db_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"  Fixed {entries_fixed} entries")
    return entries_fixed


def main():
    """Fix all rule databases."""
    print("Fixing Rule Database Reagents")
    print("=" * 70)
    print()
    
    db_dir = Path("data/conditionDB")
    
    total_fixed = 0
    
    # Fix Suzuki database
    suzuki_path = db_dir / "Suzuki_db.json"
    if suzuki_path.exists():
        total_fixed += fix_suzuki_reagents(str(suzuki_path))
        print()
    
    # Fix C-N coupling databases
    for cn_db in ["C_N_Coupling_Cu_db.json", "C_N_Coupling_Pd_db.json", "C_N_Coupling_Ni_db.json"]:
        cn_path = db_dir / cn_db
        if cn_path.exists():
            total_fixed += fix_cn_coupling_reagents(str(cn_path))
            print()
    
    print("=" * 70)
    print(f"Total entries fixed: {total_fixed}")
    print("Done!")
    
    return 0


if __name__ == "__main__":
    exit(main())
