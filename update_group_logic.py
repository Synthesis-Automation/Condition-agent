#!/usr/bin/env python3
"""Update group_logic.json to use new group IDs with '-' prefix for substituents."""

import json
from pathlib import Path

# Define the mapping from old IDs to new IDs
ID_MAPPING = {
    # Halogens
    "F": "-F", "Cl": "-Cl", "Br": "-Br", "I": "-I",
    
    # Sulfonyl halides
    "SO2F": "-SO2F", "SO2Cl": "-SO2Cl", "SO2Br": "-SO2Br", "SO2I": "-SO2I",
    
    # Acyl halides
    "COF": "-COF", "COCl": "-COCl", "COBr": "-COBr", "COI": "-COI",
    
    # Sulfonates
    "OTf": "-OTf", "OTs": "-OTs", "OMs": "-OMs", "OSO2R": "-OSO2R",
    
    # Amines
    "NH2": "-NH2", "NHR": "-NHR", "NR2": "-NR2",
    
    # Amides
    "CONH2": "-CONH2", "CONHR": "-CONHR", "CONR2": "-CONR2", "LactamN": "-LactamN",
    
    # Sulfonamides
    "SO2NH2": "-SO2NH2", "SO2NHR": "-SO2NHR", "SO2NR2": "-SO2NR2",
    
    # Alcohols/ethers
    "OH": "-OH", "OR": "-OR", "SH": "-SH", "SR": "-SR",
    
    # Carbamates
    "NCO2R": "-NCO2R",
    
    # Fluorinated
    "CF3": "-CF3", "OCF3": "-OCF3",
    
    # Carbonyls
    "CHO": "-CHO", "COR": "-COR", "CO2H": "-CO2H", "CO2R": "-CO2R",
    
    # Alkyl groups
    "CH3": "-CH3", "RCH2": "-RCH2", "R2CH": "-R2CH", "R3C": "-R3C", 
    "Bn": "-Bn", "Allyl": "-Allyl",
    
    # Organoboron
    "B(OH)2": "-B(OH)2", "Bpin": "-Bpin", "B(OR)2": "-B(OR)2", "BF3K": "-BF3K",
    
    # Organometallics (need * suffix too)
    "Sn": "-Sn*", "Zn": "-Zn*", "Mg": "-Mg*", "Si": "-Si*", "M": "-M*",
}

def update_members(members, mapping):
    """Update a members list using the ID mapping."""
    updated = []
    for member in members:
        if member in mapping:
            updated.append(mapping[member])
        else:
            updated.append(member)  # Keep as-is (scaffolds or group set names)
    return updated

def main():
    logic_file = Path("chemtools/taxonomy/data/group_logic.json")
    
    print(f"Reading {logic_file}...")
    with open(logic_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    total_updates = 0
    
    # Update all group sets
    for set_name, set_data in data.get("group_sets", {}).items():
        if "members" in set_data:
            old_members = set_data["members"]
            new_members = update_members(old_members, ID_MAPPING)
            
            if old_members != new_members:
                changes = sum(1 for old, new in zip(old_members, new_members) if old != new)
                print(f"  {set_name}: Updated {changes} members")
                print(f"    Old: {old_members}")
                print(f"    New: {new_members}")
                set_data["members"] = new_members
                total_updates += changes
    
    # Write back
    print(f"\nWriting updated {logic_file}...")
    with open(logic_file, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    
    print(f"\n✓ Updated {total_updates} group ID references in group_logic.json")

if __name__ == "__main__":
    main()
