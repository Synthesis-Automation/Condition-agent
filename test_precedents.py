#!/usr/bin/env python3
"""Test script to check precedent data in recommendations."""

import sys
import json
from chemtools import chem

# Get recommendations
result = chem.recommend.conditions(
    reaction='Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
    k=20,
    search_all_families=True
)

recommendations = result.get('recommendations', [])
print(f"Total recommendations: {len(recommendations)}\n")

for idx, rec in enumerate(recommendations, 1):
    print(f"{'='*80}")
    print(f"Recommendation #{idx}")
    print(f"{'='*80}")
    
    summary = rec.get('summary', {})
    
    # Check support
    support = summary.get('support', {})
    print(f"Support: {json.dumps(support, indent=2)}")
    
    # Check precedents
    precedents = summary.get('precedents', [])
    print(f"\nPrecedents count: {len(precedents)}")
    
    if precedents:
        print("Precedents data:")
        for i, prec in enumerate(precedents, 1):
            print(f"  {i}. {prec.get('reaction_id')}")
    else:
        print("Precedents: EMPTY LIST")
    
    print()
