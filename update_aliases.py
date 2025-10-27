#!/usr/bin/env python3
"""Update aliases.json to use new reactant_type IDs."""

import json

# Load aliases
with open('chemtools/taxonomy/data/aliases.json', 'r') as f:
    aliases = json.load(f)

# Mapping from old IDs to new IDs
mappings = {
    'alkyne': 'Alkyne',
    'acyl_source': 'RCO2H',
    'aldehyde': 'Aldehyde',
    'aliphatic_amine': 'RNH2/R2NH',
    'aniline_type': 'ArNH2/Ar2NH',
    'arx': 'ArX*',
    'alkyl_x': 'Alkyl-X',
    'roh': 'ROH',
    'aroh': 'ArOH',
    'rsh': 'RSH',
    'grignard': 'RMgX',
    'organolithium': 'RLi',
    'organozinc': 'RZnX',
    'metal_hydride': 'Metal-H',
    'hydrogen_source': 'H2',
    'diene': 'Diene',
    'dienophile': 'Dienophile',
    'oxidant': 'Oxidant',
    'borane': 'Borane',
    'ester': 'Ester',
    'enolate': 'Enolate',
    'cyanide': 'CN-',
    'iodide_nucleophile': 'I-',
    'alkoxide': 'RO-',
    'ketone': 'Ketone',
    'rb': 'RB*',
    'arb': 'ArB*',
    'vinylx': 'VinylX*',
    'heterocyclic_halide': 'HetAr-X',
    'organometallic': 'R-M',
    'alkyl_c_h': 'Alkyl-C-H',
    'arh': 'ArH',
    'amide_type': 'Amide',
    'rco2h': 'RCO2H',
    'rco2h_or_activated_acyl': 'RCO2H-activated',
    'ylide_phosphonate': 'Ylide',
    'enamine': 'Enamine',
    'imine': 'Imine',
    'azide': 'Azide',
    'nitrile': 'Nitrile',
    'epoxide': 'Epoxide',
    'peroxide': 'Peroxide',
    'carbanion': 'Carbanion',
    'alkene': 'alkene',  # This one stays lowercase
}

# Update aliases
updated_count = 0
for entry in aliases:
    if entry.get('entity_type') == 'reactant_type':
        old_id = entry['entity_id']
        if old_id in mappings:
            entry['entity_id'] = mappings[old_id]
            print(f"Updated: {old_id} -> {mappings[old_id]}")
            updated_count += 1

print(f"\nTotal updated: {updated_count}")

# Save updated aliases
with open('chemtools/taxonomy/data/aliases.json', 'w') as f:
    json.dump(aliases, f, indent=2)

print("✓ aliases.json updated successfully")
