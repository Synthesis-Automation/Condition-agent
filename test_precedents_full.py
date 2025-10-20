#!/usr/bin/env python3
"""Test script to check ALL precedent data locations."""

import sys
import json
from chemtools import chem

# Get recommendations
result = chem.recommend.conditions(
    reaction='Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
    k=20,
    search_all_families=True
)

print("=== Top-level result keys ===")
print(list(result.keys()))

print("\n=== Result 'precedents' field ===")
precedents = result.get('precedents', {})
print(f"Precedents type: {type(precedents)}")
print(f"Precedents keys: {list(precedents.keys()) if isinstance(precedents, dict) else 'N/A'}")
if isinstance(precedents, dict):
    for key, value in list(precedents.items())[:5]:
        print(f"  {key}: {value}")

print("\n=== Result 'precedents_used' field ===")
precedents_used = result.get('precedents_used', [])
print(f"Precedents used type: {type(precedents_used)}")
if isinstance(precedents_used, dict):
    print(f"Keys: {list(precedents_used.keys())}")
    print(f"Content: {precedents_used}")
elif isinstance(precedents_used, list):
    print(f"Count: {len(precedents_used)}")
    for i, prec_id in enumerate(precedents_used[:5], 1):
        print(f"  {i}. {prec_id}")
else:
    print(f"Value: {precedents_used}")

recommendations = result.get('recommendations', [])
print(f"\n=== Recommendations: {len(recommendations)} ===")

for idx, rec in enumerate(recommendations, 1):
    print(f"\nRecommendation #{idx}:")
    summary = rec.get('summary', {})
    
    # Core
    core = summary.get('core', 'N/A')
    print(f"  Core: {core}")
    
    # Support
    support = summary.get('support', {})
    count = support.get('count', 0)
    print(f"  Support count: {count}")
    
    # Precedents in summary
    precedents = summary.get('precedents', [])
    print(f"  Precedents in summary: {len(precedents)}")
    
    # Combo
    combo = rec.get('combo', {})
    print(f"  Combo: {combo}")
