#!/usr/bin/env python
"""Find and analyze specific Heck, Kumada, Negishi reactions in UNKNOWN."""

import sys
sys.path.insert(0, '.')

import importlib.util
spec = importlib.util.spec_from_file_location("sample_reactions", "tests/sample_reactions.py")
sample_reactions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sample_reactions)
SAMPLE_REACTIONS = sample_reactions.SAMPLE_REACTIONS
BUCHWALD_HARTWIG_REACTIONS = sample_reactions.BUCHWALD_HARTWIG_REACTIONS

from chemtools.analysis import analyze_reaction

def extract_smiles_from_entry(entry: str) -> tuple[str, str]:
    last_paren = entry.rfind(" (")
    if last_paren != -1:
        smiles = entry[:last_paren].strip()
        description = entry[last_paren+2:-1].strip()
        return smiles, description
    return entry.strip(), ""

all_reactions = list(SAMPLE_REACTIONS) + list(BUCHWALD_HARTWIG_REACTIONS)

heck_count = 0
kumada_count = 0
negishi_count = 0

for i, entry in enumerate(all_reactions, 1):
    smiles, desc = extract_smiles_from_entry(entry)
    if not smiles or ">>" not in smiles:
        continue
    
    desc_lower = desc.lower()
    
    if "heck" in desc_lower and heck_count < 3:
        result = analyze_reaction(smiles, use_rxn_insight=False)
        family = result.get('family', {}).get('canonical_id', 'UNKNOWN')
        hits = result.get('family', {}).get('detected', {}).get('hits', {})
        
        heck_count += 1
        print(f"\nHeck #{heck_count}")
        print(f"  Description: {desc}")
        print(f"  SMILES: {smiles}")
        print(f"  Detected: {family}")
        print(f"  Hits: aryl_halide={hits.get('aryl_halide')}, alkene={hits.get('alkene')}, boron={hits.get('boron')}")
        
    elif "kumada" in desc_lower and kumada_count < 3:
        result = analyze_reaction(smiles, use_rxn_insight=False)
        family = result.get('family', {}).get('canonical_id', 'UNKNOWN')
        hits = result.get('family', {}).get('detected', {}).get('hits', {})
        
        kumada_count += 1
        print(f"\nKumada #{kumada_count}")
        print(f"  Description: {desc}")
        print(f"  SMILES: {smiles}")
        print(f"  Detected: {family}")
        print(f"  Hits: aryl_halide={hits.get('aryl_halide')}, grignard={hits.get('grignard')}, carbonyl={hits.get('carbonyl')}")
        
    elif "negishi" in desc_lower and negishi_count < 3:
        result = analyze_reaction(smiles, use_rxn_insight=False)
        family = result.get('family', {}).get('canonical_id', 'UNKNOWN')
        hits = result.get('family', {}).get('detected', {}).get('hits', {})
        
        negishi_count += 1
        print(f"\nNegishi #{negishi_count}")
        print(f"  Description: {desc}")
        print(f"  SMILES: {smiles}")
        print(f"  Detected: {family}")
        print(f"  Hits: aryl_halide={hits.get('aryl_halide')}, organozinc={hits.get('organozinc')}, carbonyl={hits.get('carbonyl')}")

print(f"\n\nTotal found: Heck={heck_count}, Kumada={kumada_count}, Negishi={negishi_count}")
