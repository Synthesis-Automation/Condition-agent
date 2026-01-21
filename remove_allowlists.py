#!/usr/bin/env python3
"""Remove allowlists from reagent_roles.v2.json, preserving CAS numbers as notes."""

import json
from pathlib import Path

def remove_allowlists():
    file_path = Path("chemtools/taxonomy/data/reagent_roles.v2.json")
    
    print(f"Reading {file_path}...")
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    families_updated = 0
    cas_preserved = 0
    
    for family in data.get("families", []):
        if "allowlists" in family:
            al = family["allowlists"]
            
            # Check if there are CAS numbers to preserve
            cas_list = al.get("cas", [])
            if cas_list and len(cas_list) > 0:
                # Add CAS numbers to existing notes/description
                cas_str = ", ".join(cas_list)
                if "notes" in family and family["notes"]:
                    family["notes"] += f" CAS: {cas_str}"
                elif "description" in family and family["description"]:
                    family["description"] += f" (CAS: {cas_str})"
                else:
                    family["notes"] = f"CAS: {cas_str}"
                
                print(f"  {family['id']}: Preserved {len(cas_list)} CAS numbers")
                cas_preserved += len(cas_list)
            
            # Remove allowlists
            del family["allowlists"]
            families_updated += 1
    
    # Write back
    print(f"\nWriting updated {file_path}...")
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    
    print(f"\n✓ Removed allowlists from {families_updated} families")
    print(f"✓ Preserved {cas_preserved} CAS numbers in notes/description")

if __name__ == "__main__":
    remove_allowlists()
