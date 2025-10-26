#!/usr/bin/env python
"""Analyze specific UNKNOWN reactions to see why they're not detected."""

import sys
sys.path.insert(0, '.')

from chemtools.analysis import analyze_reaction
import re

# Read sample_reactions.py to extract specific reactions
with open('tests/sample_reactions.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Parse entries
entries = re.findall(r'\((r?"[^"]+",\s*"[^"]+",\s*"[^"]+")\)', content)

unknown_count = 0
heck_found = 0
kumada_found = 0
negishi_found = 0

for entry in entries:
    parts = entry.strip('()').split(',')
    smiles = parts[0].strip().strip('"').strip("r'").strip("'")
    
    result = analyze_reaction(smiles, use_rxn_insight=False)
    detection = result.get('family', {}).get('detected', {})
    family = detection.get('family', 'UNKNOWN')
    
    if family == 'UNKNOWN':
        unknown_count += 1
        desc = parts[1].strip().strip('"')
        
        # Look for Heck/Kumada/Negishi
        if 'heck' in desc.lower():
            heck_found += 1
            if heck_found <= 2:  # Show first 2
                hits = detection.get('hits', {})
                print(f"\nHeck reaction #{heck_found}:")
                print(f"  Description: {desc}")
                print(f"  SMILES: {smiles}")
                print(f"  Detected: {family}")
                print(f"  Key hits:")
                for k in ['aryl_halide', 'alkene', 'boron', 'grignard']:
                    print(f"    {k}: {hits.get(k, False)}")
                    
        elif 'kumada' in desc.lower():
            kumada_found += 1
            if kumada_found <= 2:
                hits = detection.get('hits', {})
                print(f"\nKumada reaction #{kumada_found}:")
                print(f"  Description: {desc}")
                print(f"  SMILES: {smiles}")
                print(f"  Detected: {family}")
                print(f"  Key hits:")
                for k in ['aryl_halide', 'grignard', 'carbonyl']:
                    print(f"    {k}: {hits.get(k, False)}")
                    
        elif 'negishi' in desc.lower():
            negishi_found += 1
            if negishi_found <= 2:
                hits = detection.get('hits', {})
                print(f"\nNegishi reaction #{negishi_found}:")
                print(f"  Description: {desc}")
                print(f"  SMILES: {smiles}")
                print(f"  Detected: {family}")
                print(f"  Key hits:")
                for k in ['aryl_halide', 'organozinc', 'carbonyl']:
                    print(f"    {k}: {hits.get(k, False)}")

print(f"\n\nSummary:")
print(f"  Total UNKNOWN: {unknown_count}")
print(f"  Heck in UNKNOWN: {heck_found}")
print(f"  Kumada in UNKNOWN: {kumada_found}")
print(f"  Negishi in UNKNOWN: {negishi_found}")
