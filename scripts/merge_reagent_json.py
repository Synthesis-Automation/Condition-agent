#!/usr/bin/env python3
"""
Merge reagents_export.json into reagents.csv.
- For existing CAS: keep CSV name/abbreviation/roles, overwrite smiles, add new columns
- For new CAS: add full row with empty role columns
"""

import csv
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REAGENTS_CSV = PROJECT_ROOT / "data" / "reagent_db" / "reagents.csv"
REAGENTS_JSON = PROJECT_ROOT / "data" / "reagent_db" / "reagents_export.json"
OUTPUT_CSV = PROJECT_ROOT / "data" / "reagent_db" / "reagents_merged.csv"

# New column order
NEW_COLUMNS = [
    "name", "abbreviation", "cas", "smiles", "formula", "type", 
    "density", "mw", "bp", "mp", "volatile", "viscose",
    "role_1", "family_1", "tag_1", "role_2", "family_2", "tag_2"
]

def main():
    print(f"Loading {REAGENTS_JSON}...")
    with open(REAGENTS_JSON, "r", encoding="utf-8") as f:
        json_data = json.load(f)
    
    # Index JSON by CAS
    json_by_cas = {}
    for item in json_data["data"]:
        cas = (item.get("casno") or "").strip()
        if cas:
            json_by_cas[cas] = item
    
    print(f"Loaded {len(json_by_cas)} entries from JSON")
    
    print(f"Loading {REAGENTS_CSV}...")
    csv_rows = []
    csv_cas_set = set()
    
    with open(REAGENTS_CSV, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            csv_rows.append(row)
            cas = (row.get("cas") or "").strip()
            if cas:
                csv_cas_set.add(cas)
    
    print(f"Loaded {len(csv_rows)} entries from CSV")
    
    # Merge existing entries
    merged_rows = []
    updated_count = 0
    
    for csv_row in csv_rows:
        cas = (csv_row.get("cas") or "").strip()
        
        # Start with CSV values (keep name, abbreviation, roles)
        merged = {
            "name": csv_row.get("name", ""),
            "abbreviation": csv_row.get("abbreviation", ""),
            "cas": cas,
            "role_1": csv_row.get("role_1", ""),
            "family_1": csv_row.get("family_1", ""),
            "tag_1": csv_row.get("tag_1", ""),
            "role_2": csv_row.get("role_2", ""),
            "family_2": csv_row.get("family_2", ""),
            "tag_2": csv_row.get("tag_2", ""),
        }
        
        # If CAS exists in JSON, update/add fields from JSON
        if cas in json_by_cas:
            json_item = json_by_cas[cas]
            merged["smiles"] = json_item.get("smile", "")
            merged["formula"] = json_item.get("formula", "")
            merged["type"] = json_item.get("type", "")
            merged["density"] = str(json_item.get("density", "")) if json_item.get("density") is not None else ""
            merged["mw"] = str(json_item.get("mw", "")) if json_item.get("mw") is not None else ""
            merged["bp"] = str(json_item.get("bp", "")) if json_item.get("bp") is not None else ""
            merged["mp"] = str(json_item.get("mp", "")) if json_item.get("mp") is not None else ""
            merged["volatile"] = str(json_item.get("volatile", "")) if json_item.get("volatile") is not None else ""
            merged["viscose"] = str(json_item.get("viscose", "")) if json_item.get("viscose") is not None else ""
            updated_count += 1
        else:
            # No JSON match - keep existing smiles or empty for new columns
            merged["smiles"] = csv_row.get("smiles", "")
            merged["formula"] = ""
            merged["type"] = ""
            merged["density"] = ""
            merged["mw"] = ""
            merged["bp"] = ""
            merged["mp"] = ""
            merged["volatile"] = ""
            merged["viscose"] = ""
        
        merged_rows.append(merged)
    
    print(f"Updated {updated_count} existing entries")
    
    # Add new entries from JSON
    new_count = 0
    for cas, json_item in json_by_cas.items():
        if cas not in csv_cas_set:
            new_row = {
                "name": json_item.get("name", ""),
                "abbreviation": json_item.get("abbreviation", ""),
                "cas": cas,
                "smiles": json_item.get("smile", ""),
                "formula": json_item.get("formula", ""),
                "type": json_item.get("type", ""),
                "density": str(json_item.get("density", "")) if json_item.get("density") is not None else "",
                "mw": str(json_item.get("mw", "")) if json_item.get("mw") is not None else "",
                "bp": str(json_item.get("bp", "")) if json_item.get("bp") is not None else "",
                "mp": str(json_item.get("mp", "")) if json_item.get("mp") is not None else "",
                "volatile": str(json_item.get("volatile", "")) if json_item.get("volatile") is not None else "",
                "viscose": str(json_item.get("viscose", "")) if json_item.get("viscose") is not None else "",
                "role_1": "",
                "family_1": "",
                "tag_1": "",
                "role_2": "",
                "family_2": "",
                "tag_2": "",
            }
            merged_rows.append(new_row)
            new_count += 1
    
    print(f"Added {new_count} new entries from JSON")
    
    # Write merged CSV
    print(f"Writing to {OUTPUT_CSV}...")
    with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=NEW_COLUMNS)
        writer.writeheader()
        writer.writerows(merged_rows)
    
    print(f"Done! Total entries: {len(merged_rows)}")
    print(f"  - Updated: {updated_count}")
    print(f"  - Added: {new_count}")
    print(f"  - Unchanged: {len(csv_rows) - updated_count}")
    print(f"\nMerged file saved to: {OUTPUT_CSV}")
    print("Review the file, then replace reagents.csv if satisfied.")

if __name__ == "__main__":
    main()
