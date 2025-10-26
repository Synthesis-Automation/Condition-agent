#!/usr/bin/env python
"""Test hydrogenation detection after fix."""

from chemtools.analysis import analyze_reaction

r = analyze_reaction('c1ccc(C=O)cc1.[H][H]>>c1ccc(CO)cc1', use_rxn_insight=False)
family = r.get('family', {}).get('canonical_id', 'UNKNOWN')
confidence = r.get('family', {}).get('detected', {}).get('confidence', 0.0)
hits = r.get('family', {}).get('detected', {}).get('hits', {})

print(f"Family: {family}")
print(f"Confidence: {confidence}")
print(f"Carbonyl hit: {hits.get('carbonyl', False)}")
print(f"Alkene hit: {hits.get('alkene', False)}")
