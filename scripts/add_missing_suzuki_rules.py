"""Add 8 missing Suzuki database entries that were lost during git checkout."""

import json

# Load current database
with open('data/conditionDB/suzuki_db.json', 'r', encoding='utf-8') as f:
    db = json.load(f)

print(f"Current entries: {len(db['entries'])}")

# Define the 8 missing entries
new_entries = [
    {
        "id": "SCDB-SUZ-BULKY-NUC-XPHOS",
        "reaction_type": "Suzuki_Miyaura",
        "name": "Hindered nucleophile (2,6-disubstituted aryl boron) + ArBr (XPhos)",
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "token_signature": ["Suzuki", "hindered_aryl_boron", "XPhos", "Pd2(dba)3", "K3PO4"],
        "conditions": {
            "pd_source": ["Pd2(dba)3"],
            "ligands": ["XPhos"],
            "boron_partner": ["2,6-disubstituted aryl-B(OH)2 1.3 eq"],
            "base": ["K3PO4 2.5 eq"],
            "solvent": ["toluene"],
            "temperature_C": [80, 90],
            "time_h": [12, 18],
            "loadings": {"Pd_mol%": [2.0, 3.0], "ligand_mol%": [4.0, 6.0], "base_eq": [2.5, 3.0]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "Br", "nucleophile.ortho_sub_count": "2"},
            "feature_requirements": {"electrophile.lg_class": ["Br"], "nucleophile.ortho_sub_count": [2]}
        },
        "evidence": {"provenance": "literature_theme", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["Bulky XPhos ligand handles sterically hindered aryl boron coupling partners."],
        "priority": 48
    },
    {
        "id": "SCDB-SUZ-2HETARYL-SPHOS",
        "reaction_type": "Suzuki_Miyaura",
        "name": "2-Heteroaryl halides (N-adjacent position) + aryl boron (SPhos)",
        "rxn_smiles_min": "[c,n:1]-[Br,I:2].[B:3](-[O])(-[O])-[c:4]>>[c,n:1]-[c:4]",
        "token_signature": ["Suzuki", "2-hetaryl_halide", "SPhos", "Pd(OAc)2", "K2CO3"],
        "conditions": {
            "pd_source": ["Pd(OAc)2"],
            "ligands": ["SPhos"],
            "boron_partner": ["aryl-B(OH)2 1.3 eq"],
            "base": ["K2CO3 2.0 eq"],
            "solvent": ["THF/H2O (9:1)"],
            "temperature_C": [50, 65],
            "time_h": [8, 14],
            "loadings": {"Pd_mol%": [3.0, 5.0], "ligand_mol%": [6.0, 10.0], "base_eq": [2.0, 2.5]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "I/Br", "electrophile.ring_hetero_count": "≥1"},
            "feature_requirements": {"electrophile.lg_class": ["I", "Br"], "electrophile.ring_hetero_count": [1, 2, 3]}
        },
        "evidence": {"provenance": "literature_theme", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["2-Heteroaryl positions prone to reduced reactivity; SPhos helps."],
        "priority": 65
    },
    {
        "id": "SCDB-SUZ-ARI-META-TBUXPHOS",
        "reaction_type": "Suzuki_Miyaura",
        "name": "meta-Hindered ArI (3,5-disubstituted) + aryl boron (tBuXPhos)",
        "rxn_smiles_min": "[c:1]-[I:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "token_signature": ["Suzuki", "meta_hindered_ArI", "tBuXPhos", "Pd(OAc)2", "K3PO4"],
        "conditions": {
            "pd_source": ["Pd(OAc)2"],
            "ligands": ["tBuXPhos"],
            "boron_partner": ["aryl-B(OH)2 1.3 eq"],
            "base": ["K3PO4 2.5 eq"],
            "solvent": ["THF"],
            "temperature_C": [70, 80],
            "time_h": [10, 16],
            "loadings": {"Pd_mol%": [2.0, 4.0], "ligand_mol%": [4.0, 8.0], "base_eq": [2.5, 3.0]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "I", "electrophile.meta_sub_count": "≥2"},
            "feature_requirements": {"electrophile.lg_class": ["I"], "electrophile.meta_sub_count": [2]}
        },
        "evidence": {"provenance": "literature_theme", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["Bulky tBuXPhos for meta-substitution patterns."],
        "priority": 55
    },
    {
        "id": "M-SUZ-OTf-DPPF",
        "reaction_type": "Suzuki_Miyaura",
        "name": "Aryl triflate + aryl boron (scheme-based, dppf ligand)",
        "rxn_smiles_min": "[c:1]-[#16:2](=O)(=O)-[#8]-[#6](F)(F)F.[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "token_signature": ["Suzuki", "ArOTf", "scheme_match", "PdCl2(dppf)", "K2CO3"],
        "conditions": {
            "pd_source": ["PdCl2(dppf)"],
            "ligands": ["dppf"],
            "boron_partner": ["aryl-B(OH)2 1.3 eq"],
            "base": ["K2CO3 2.0 eq"],
            "solvent": ["THF"],
            "temperature_C": [60, 70],
            "time_h": [8, 12],
            "loadings": {"Pd_mol%": [3.0, 5.0], "ligand_mol%": [0.0, 0.0], "base_eq": [2.0, 2.5]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "OTf"},
            "feature_requirements": {"electrophile.lg_class": ["OTf"]}
        },
        "evidence": {"provenance": "scheme_match", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["Scheme-based match for OTf with dppf."],
        "priority": 80
    },
    {
        "id": "M-SUZ-VINYL-RT",
        "reaction_type": "Suzuki_Miyaura",
        "name": "Vinyl halide + vinyl boron (scheme-based, RT)",
        "rxn_smiles_min": "[C:1]=[C:2]-[I,Br:3].[B:4](-[O])(-[O])-[C:5]=[C:6]>>[C:1]=[C:2]-[C:5]=[C:6]",
        "token_signature": ["Suzuki", "vinyl_coupling", "scheme_match", "Pd(PPh3)4", "RT"],
        "conditions": {
            "pd_source": ["Pd(PPh3)4"],
            "ligands": ["PPh3 (intrinsic)"],
            "boron_partner": ["vinyl-B(OH)2 1.2 eq"],
            "base": ["K2CO3 2.0 eq"],
            "solvent": ["THF/H2O (9:1)"],
            "temperature_C": [20, 30],
            "time_h": [4, 8],
            "loadings": {"Pd_mol%": [2.0, 3.0], "ligand_mol%": [0.0, 0.0], "base_eq": [2.0, 2.5]}
        },
        "env": {
            "features_from_smiles": {"electrophile.is_vinyl": "true", "nucleophile.is_vinyl": "true"},
            "feature_requirements": {"electrophile.is_vinyl": [True], "nucleophile.is_vinyl": [True]}
        },
        "evidence": {"provenance": "scheme_match", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["Vinyl-vinyl coupling at room temperature."],
        "priority": 75
    },
    {
        "id": "M-SUZ-BF3K-GENERAL",
        "reaction_type": "Suzuki_Miyaura",
        "name": "Aryl halide + BF3K boron partner (scheme-based)",
        "rxn_smiles_min": "[c:1]-[Br,I:2].[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "token_signature": ["Suzuki", "BF3K_boron", "scheme_match", "Pd(OAc)2", "SPhos"],
        "conditions": {
            "pd_source": ["Pd(OAc)2"],
            "ligands": ["SPhos"],
            "boron_partner": ["aryl-BF3K 1.2 eq"],
            "base": ["Cs2CO3 2.0 eq"],
            "solvent": ["THF/H2O (9:1)"],
            "temperature_C": [50, 65],
            "time_h": [8, 12],
            "loadings": {"Pd_mol%": [3.0, 5.0], "ligand_mol%": [6.0, 10.0], "base_eq": [2.0, 2.5]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "Br/I", "boron.partner_type": "BF3K"},
            "feature_requirements": {"electrophile.lg_class": ["Br", "I"], "boron.partner_type": ["BF3K"]}
        },
        "evidence": {"provenance": "scheme_match", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["BF3K partners are moisture-stable alternatives to boronic acids."],
        "priority": 65
    },
    {
        "id": "SCDB-SUZ-HET-AZINE-BORON",
        "reaction_type": "Suzuki_Miyaura",
        "name": "ArBr + azine-boronic acid (pyridine/pyrimidine boron partner)",
        "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])(-[O])-[c,n:4]>>[c:1]-[c,n:4]",
        "token_signature": ["Suzuki", "azine_boron_partner", "Pd(PPh3)4", "K2CO3"],
        "conditions": {
            "pd_source": ["Pd(PPh3)4"],
            "ligands": ["PPh3 (intrinsic)"],
            "boron_partner": ["azine-B(OH)2 1.3 eq"],
            "base": ["K2CO3 2.5 eq"],
            "solvent": ["dioxane/H2O (4:1)"],
            "temperature_C": [70, 85],
            "time_h": [10, 16],
            "loadings": {"Pd_mol%": [3.0, 5.0], "ligand_mol%": [0.0, 0.0], "base_eq": [2.5, 3.0]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "Br", "nucleophile.ring_hetero_count": "≥1"},
            "feature_requirements": {"electrophile.lg_class": ["Br"], "nucleophile.ring_hetero_count": [1, 2]}
        },
        "evidence": {"provenance": "literature_theme", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["Azine boron partners (pyridine, pyrimidine) require careful conditions."],
        "priority": 64
    },
    {
        "id": "SCDB-SUZ-DEFAULT-OTf",
        "reaction_type": "Suzuki_Miyaura",
        "role": "default_condition",
        "name": "General aryl triflate + aryl boron (default)",
        "rxn_smiles_min": "[c:1]-[#16:2](=O)(=O)-[#8]-[#6](F)(F)F.[B:3](-[O])(-[O])-[c:4]>>[c:1]-[c:4]",
        "token_signature": ["Suzuki", "ArOTf", "default", "PdCl2(dppf)", "K2CO3"],
        "conditions": {
            "pd_source": ["PdCl2(dppf)"],
            "ligands": ["dppf"],
            "boron_partner": ["aryl-B(OH)2 1.3 eq"],
            "base": ["K2CO3 2.0 eq"],
            "solvent": ["THF"],
            "temperature_C": [60, 70],
            "time_h": [8, 12],
            "loadings": {"Pd_mol%": [3.0, 5.0], "ligand_mol%": [0.0, 0.0], "base_eq": [2.0, 2.5]}
        },
        "env": {
            "features_from_smiles": {"electrophile.lg_class": "OTf"},
            "feature_requirements": {"electrophile.lg_class": ["OTf"]}
        },
        "evidence": {"provenance": "playbook_default", "last_updated": "2025-10-05T00:00:00Z"},
        "notes": ["Default for OTf electrophiles."],
        "priority": 0
    }
]

# Find the index of the first DEFAULT entry
default_idx = None
for i, entry in enumerate(db['entries']):
    if entry.get('role') == 'default_condition':
        default_idx = i
        break

if default_idx is None:
    print("ERROR: Could not find default entries!")
    exit(1)

print(f"Inserting 8 new entries before index {default_idx} (before {db['entries'][default_idx]['id']})")

# Insert all new entries before the defaults
for i, new_entry in enumerate(new_entries):
    db['entries'].insert(default_idx + i, new_entry)
    print(f"  Added: {new_entry['id']} (priority {new_entry.get('priority', 0)})")

print(f"\nNew total entries: {len(db['entries'])}")

# Save the updated database
with open('data/conditionDB/suzuki_db.json', 'w', encoding='utf-8') as f:
    json.dump(db, f, indent=2, ensure_ascii=False)

print("Database updated successfully!")
