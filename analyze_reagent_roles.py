#!/usr/bin/env python3
"""Analyze reagent_roles.v2.json structure."""

import json
from pathlib import Path

def analyze_reagent_roles():
    file_path = Path("chemtools/taxonomy/data/reagent_roles.v2.json")
    
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    roles = data.get("roles", [])
    families = data.get("families", [])
    
    print(f"Total roles: {len(roles)}")
    print(f"Total families: {len(families)}")
    
    # Check roles
    roles_with_name = 0
    roles_name_equals_id = 0
    for role in roles:
        if "name" in role:
            roles_with_name += 1
            # Check if name is just title case of id
            expected_name = role["id"].replace("_", " ").title()
            if role["name"] == expected_name:
                roles_name_equals_id += 1
            else:
                print(f"\nRole ID: {role['id']}")
                print(f"  Name: {role['name']}")
                print(f"  Expected: {expected_name}")
    
    print(f"\nRoles with name field: {roles_with_name}/{len(roles)}")
    print(f"Roles where name = title_case(id): {roles_name_equals_id}/{roles_with_name}")
    
    # Check families
    families_with_name = 0
    families_name_equals_id = 0
    different_names = []
    
    for family in families:
        if "name" in family:
            families_with_name += 1
            # Check if name equals id (families use snake_case for both)
            if family["name"] == family["id"]:
                families_name_equals_id += 1
            elif family["name"].replace(" ", "_") == family["id"]:
                families_name_equals_id += 1
            else:
                different_names.append((family["id"], family["name"]))
    
    print(f"\nFamilies with name field: {families_with_name}/{len(families)}")
    print(f"Families where name matches id: {families_name_equals_id}/{families_with_name}")
    print(f"Families with different names: {len(different_names)}")
    
    if different_names:
        print("\nAll families with different names:")
        for fid, fname in different_names:
            print(f"  {fid} → {fname}")

if __name__ == "__main__":
    analyze_reagent_roles()
