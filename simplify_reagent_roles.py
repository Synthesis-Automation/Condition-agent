#!/usr/bin/env python3
"""Simplify reagent_roles.v2.json by removing redundant name fields."""

import json
from pathlib import Path

def simplify_reagent_roles():
    file_path = Path("chemtools/taxonomy/data/reagent_roles.v2.json")
    
    print(f"Reading {file_path}...")
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    roles_removed = 0
    families_removed = 0
    families_aliased = 0
    
    # Process roles - remove all name fields (all are title_case of id)
    for role in data.get("roles", []):
        if "name" in role:
            print(f"  Role: Removing name '{role['name']}' from {role['id']}")
            del role["name"]
            roles_removed += 1
    
    # Process families
    for family in data.get("families", []):
        if "name" in family:
            fid = family["id"]
            fname = family["name"]
            
            # Check if name matches id (with space/underscore conversion)
            if fname == fid or fname.replace(" ", "_") == fid:
                # Redundant name - just remove it
                del family["name"]
                families_removed += 1
            else:
                # Different name - preserve as alias
                print(f"  Family: Converting name to alias for {fid}: '{fname}'")
                
                # Add to aliases if not already there
                if "aliases" not in family:
                    family["aliases"] = []
                if fname not in family["aliases"]:
                    family["aliases"].insert(0, fname)  # Put at front
                
                # Remove name field
                del family["name"]
                families_aliased += 1
    
    # Write back
    print(f"\nWriting updated {file_path}...")
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    
    print(f"\n✓ Simplified reagent_roles.v2.json:")
    print(f"  - Removed {roles_removed} role name fields")
    print(f"  - Removed {families_removed} redundant family name fields")
    print(f"  - Converted {families_aliased} family names to aliases")
    print(f"  - Total: {roles_removed + families_removed + families_aliased} name fields removed")

if __name__ == "__main__":
    simplify_reagent_roles()
