"""Fix SMARTS patterns that are matching too broadly."""

import json

# Load database
with open('data/conditionDB/suzuki_db.json', 'r', encoding='utf-8') as f:
    db = json.load(f)

# Fixes needed:
fixes = {
    "SCDB-SUZ-BULKY-NUC-XPHOS": {
        # Should only match when BOTH ortho positions on nucleophile (boron partner) are substituted
        # AND electrophile is simple ArBr (no heteroatoms, no ortho subs on electrophile)
        "rxn_smiles_min": "[c:1]1[cH][cH][cH][cH][cH]1-[Br:2].[B:3](-[O])(-[O])-[c:4]([#6,#7,#8,#16])([#6,#7,#8,#16])>>[c:1]-[c:4]",
        "reason": "Electrophile must be SIMPLE phenyl (all [cH]), nucleophile has both ortho subs"
    },
    "SCDB-SUZ-ARBR-ORTHO-XPhos": {
        # Should match when electrophile has ortho substituent
        "rxn_smiles_min": "[c:1]1[c]([#6,#7,#8,#16])[cH][cH][cH][cH]1-[Br:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "reason": "Electrophile has ortho sub, nucleophile is simple"
    },
    "SCDB-SUZ-2HETARYL-SPHOS": {
        # Should match heteroaromatic electrophile (with N in ring)
        "rxn_smiles_min": "[c:1]1[n,o,s][c,n][c,n][c,n][c,n]1-[Br,I:2].[B:3](-[O])(-[O])-[c:4]1[cH][cH][cH][cH][cH]1>>[c:1]-[c:4]",
        "reason": "Electrophile is heteroaromatic, nucleophile is simple phenyl"
    },
    "SCDB-SUZ-HET-2PYRIDYL-SLOWRELEASE": {
        # BF3K pattern - use simpler approach (just match any ArBr + pyridyl boron)
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3]-[c:4]1[n][c,n][c,n][c,n][c,n]1>>[c:1]-[c:4]",
        "reason": "Match ArBr + pyridyl-boron (BF3K would be hard to represent in product SMARTS)"
    },
    "SCDB-SUZ-VINYL-DPPF-RT": {
        # Keep vinyl pattern but increase priority to beat M-SUZ-VINYL-RT if needed
        # Actually, these are duplicates - let's keep M-SUZ-VINYL-RT (priority 75) as the main one
        # and make this one more specific OR accept it won't be matched
        "priority": 72,  # Still lower than M-SUZ (75), so this will be beaten - that's OK
        "reason": "Lower priority than M-SUZ-VINYL-RT - expected to be overridden"
    }
}

for entry in db['entries']:
    entry_id = entry['id']
    
    if entry_id in fixes:
        fix = fixes[entry_id]
        
        print(f"Fixing {entry_id}:")
        
        if 'rxn_smiles_min' in fix:
            print(f"  Old SMARTS: {entry['rxn_smiles_min']}")
            print(f"  New SMARTS: {fix['rxn_smiles_min']}")
            entry['rxn_smiles_min'] = fix['rxn_smiles_min']
        
        if 'priority' in fix:
            print(f"  Old priority: {entry.get('priority', 'N/A')}")
            print(f"  New priority: {fix['priority']}")
            entry['priority'] = fix['priority']
        
        print(f"  Reason: {fix['reason']}")
        print()

# Save
with open('data/conditionDB/suzuki_db.json', 'w', encoding='utf-8') as f:
    json.dump(db, f, indent=2, ensure_ascii=False)

print("Database updated!")
