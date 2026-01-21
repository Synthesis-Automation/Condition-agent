#!/usr/bin/env python3
"""Analyze allowlists usage in reagent_roles.v2.json."""

import json
from pathlib import Path

def analyze_allowlists():
    file_path = Path("chemtools/taxonomy/data/reagent_roles.v2.json")
    
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    families = data.get("families", [])
    
    total = len(families)
    with_allowlists = 0
    non_empty_cas = 0
    non_empty_names = 0
    non_empty_keywords = 0
    
    cas_entries = []
    name_entries = []
    
    for family in families:
        if "allowlists" in family:
            with_allowlists += 1
            al = family["allowlists"]
            
            if al.get("cas") and len(al["cas"]) > 0:
                non_empty_cas += 1
                cas_entries.append((family["id"], len(al["cas"])))
            
            if al.get("names") and len(al["names"]) > 0:
                non_empty_names += 1
                name_entries.append((family["id"], len(al["names"])))
            
            if al.get("keywords") and len(al["keywords"]) > 0:
                non_empty_keywords += 1
    
    print(f"Total families: {total}")
    print(f"Families with allowlists: {with_allowlists}")
    print(f"\nNon-empty arrays:")
    print(f"  CAS numbers: {non_empty_cas}")
    print(f"  Names: {non_empty_names}")
    print(f"  Keywords: {non_empty_keywords}")
    
    if cas_entries:
        print(f"\nFamilies with CAS numbers ({len(cas_entries)}):")
        for fid, count in cas_entries:
            print(f"  {fid}: {count} CAS numbers")
    
    if name_entries:
        print(f"\nFamilies with names ({len(name_entries)}):")
        for fid, count in name_entries:
            print(f"  {fid}: {count} names")

if __name__ == "__main__":
    analyze_allowlists()
