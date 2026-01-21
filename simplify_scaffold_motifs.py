#!/usr/bin/env python3
"""Remove redundant 'name' field from scaffold_motifs.v1.3.json."""

import json
from pathlib import Path

def main():
    file_path = Path("chemtools/taxonomy/data/scaffold_motifs.v1.3.json")
    
    print(f"Reading {file_path}...")
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    removed_count = 0
    
    # Process each compound
    for entry in data.get("compounds", []):
        if "name" in entry and "id" in entry:
            if entry["name"] == entry["id"]:
                print(f"  Removing redundant name from: {entry['id']}")
                del entry["name"]
                removed_count += 1
            else:
                # Name differs from ID - preserve as alias
                name_value = entry["name"]
                print(f"  Converting name to alias for {entry['id']}: '{name_value}'")
                
                # Add to aliases if not already there
                if "aliases" not in entry:
                    entry["aliases"] = []
                if name_value not in entry["aliases"]:
                    entry["aliases"].append(name_value)
                
                # Remove name field
                del entry["name"]
                removed_count += 1
    
    # Write back
    print(f"\nWriting updated {file_path}...")
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    
    print(f"\n✓ Removed {removed_count} redundant 'name' fields from scaffold motifs")

if __name__ == "__main__":
    main()
