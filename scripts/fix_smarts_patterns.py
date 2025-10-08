"""Update SMARTS patterns to be more specific and remove unnecessary feature requirements."""

import json

# Load database
with open('data/conditionDB/suzuki_db.json', 'r', encoding='utf-8') as f:
    db = json.load(f)

print(f"Current entries: {len(db['entries'])}\n")

# Define SMARTS pattern updates
smarts_updates = {
    "SCDB-SUZ-BULKY-NUC-XPHOS": {
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]([#6,#7,#8,#16])([#6,#7,#8,#16])>>[c:1]-[c:4]",
        "remove_features": True,
        "reason": "SMARTS now requires both ortho positions substituted (not H)"
    },
    "SCDB-SUZ-2HETARYL-SPHOS": {
        "rxn_smiles_min": "[c:1]1[n,o,s][c,n][c,n][c,n][c,n]1-[Br,I:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "remove_features": True,
        "reason": "SMARTS now explicitly matches heteroaromatic rings (5-6 membered with N/O/S)"
    },
    "SCDB-SUZ-ARI-META-TBUXPHOS": {
        "rxn_smiles_min": "[c:1]1[c]([#6,#7,#8,#16])[c][c]([#6,#7,#8,#16])[c][c]1-[I:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "remove_features": True,
        "reason": "SMARTS now requires meta positions (3,5) substituted"
    },
    "SCDB-SUZ-ARBR-ORTHO-XPhos": {
        "rxn_smiles_min": "[c:1]1[c]([#6,#7,#8,#16])[c][c][c][c]1-[Br:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "remove_features": True,
        "reason": "SMARTS now requires ortho position substituted"
    },
    "SCDB-SUZ-HET-2PYRIDYL-SLOWRELEASE": {
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[F])(-[F])(-[F])-[K].[c:4]1[n][c,n][c,n][c,n][c,n]1>>[c:1]-[c:4]",
        "remove_features": False,
        "reason": "SMARTS for BF3K salt form (slow-release boron)"
    },
    "SCDB-SUZ-HET-AZINE-BORON": {
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]1[n][c,n][c,n][c,n][c,n]1>>[c:1]-[c:4]",
        "remove_features": True,
        "reason": "SMARTS explicitly matches azine (pyridine/pyrimidine) on boron side"
    },
    "SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE": {
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]1[o,s][c][c][c]1>>[c:1]-[c:4]",
        "remove_features": True,
        "reason": "SMARTS explicitly matches thiophene/furan on boron side"
    },
    "SCDB-SUZ-VINYL-DPPF-RT": {
        "rxn_smiles_min": "[C:1]=[C:2]-[Br:3].[B:4](-[O])(-[O])-[C:5]=[C:6]>>[C:1]=[C:2]-[C:5]=[C:6]",
        "remove_features": True,
        "reason": "SMARTS already specific for vinyl-vinyl coupling"
    },
    "SCDB-SUZ-ALKYL-PRIMARYI-9BBN": {
        "rxn_smiles_min": "[C;X4;H2:1]-[I:2].[B:3]>>[C:1]-[B:3]",
        "remove_features": True,
        "reason": "SMARTS now requires primary sp3 carbon (X4=4 bonds, H2=2 hydrogens)"
    }
}

# Apply updates
updated_count = 0
for entry in db['entries']:
    entry_id = entry['id']
    
    if entry_id in smarts_updates:
        update = smarts_updates[entry_id]
        
        print(f"Updating {entry_id}:")
        print(f"  Old SMARTS: {entry['rxn_smiles_min']}")
        print(f"  New SMARTS: {update['rxn_smiles_min']}")
        print(f"  Reason: {update['reason']}")
        
        # Update SMARTS
        entry['rxn_smiles_min'] = update['rxn_smiles_min']
        
        # Remove feature requirements if specified
        if update['remove_features']:
            if 'env' in entry and 'feature_requirements' in entry['env']:
                old_reqs = entry['env']['feature_requirements']
                print(f"  Removed feature_requirements: {old_reqs}")
                entry['env']['feature_requirements'] = {}
        
        print()
        updated_count += 1

print(f"\nUpdated {updated_count} entries")

# Save
with open('data/conditionDB/suzuki_db.json', 'w', encoding='utf-8') as f:
    json.dump(db, f, indent=2, ensure_ascii=False)

print("Database saved!")
