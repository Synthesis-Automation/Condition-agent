#!/usr/bin/env python3
"""
Fix Suzuki Database
===================
Apply all fixes to suzuki_db.json:
1. Add priorities to entries without them
2. Fix boronic acid SMARTS patterns
"""

import json
from pathlib import Path

# Load database
db_path = Path("C:/Git-softwares/Condition-agent/data/conditionDB/suzuki_db.json")
with open(db_path, 'r', encoding='utf-8') as f:
    db = json.load(f)

# Priority assignments
priorities = {
    "SCDB-SUZ-ARBRI-GENERAL-PPh3": 45,
    "SCDB-SUZ-ARBRI-GENERAL-SPhos": 46,
    "SCDB-SUZ-ARCL-EPoor-XPhos": 60,
    "SCDB-SUZ-ARCL-ERich-L95": 58,
    "SCDB-SUZ-OTf-DPPF": 70,
    "SCDB-SUZ-HET-2PYRIDYL-SLOWRELEASE": 68,
    "SCDB-SUZ-VINYL-DPPF-RT": 72,
    "SCDB-SUZ-ALKYL-PRIMARYI-9BBN": 62,
    "SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE": 66,
    "SCDB-SUZ-MIYAURA-BORYLATION": 52,
    "SCDB-SUZ-HET-AZINE-BORON": 64,
}

# Apply fixes
for entry in db["entries"]:
    entry_id = entry.get("id")
    
    # Add priorities
    if entry_id in priorities and "priority" not in entry:
        entry["priority"] = priorities[entry_id]
        print(f"Added priority {priorities[entry_id]} to {entry_id}")
    
    # Fix boronic acid SMARTS patterns
    if "rxn_smiles_min" in entry:
        old = entry["rxn_smiles_min"]
        # Fix [B:N]([O])[O] -> [B:N](-[O])(-[O])
        new = old.replace("[B:3]([O])[O]", "[B:3](-[O])(-[O])")
        new = new.replace("[B:4]([O])[O]", "[B:4](-[O])(-[O])")
        if old != new:
            entry["rxn_smiles_min"] = new
            print(f"Fixed SMARTS in {entry_id}: {old} -> {new}")

# Save fixed database
with open(db_path, 'w', encoding='utf-8') as f:
    json.dump(db, f, indent=2, ensure_ascii=False)

print("\nDatabase fixed and saved!")
print(f"Total entries: {len(db['entries'])}")
