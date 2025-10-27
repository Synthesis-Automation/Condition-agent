#!/usr/bin/env python3
"""
Debug script to trace why detection isn't working in cross-family mode.
"""

import sys
sys.path.insert(0, '.')

# Patch the recommender to add debug output
import chemtools.recommend.modules.recommender as rec_module

original_recommend = rec_module.recommend_from_reaction

def debug_recommend(*args, **kwargs):
    """Wrapper to add debug output."""
    print("\n=== DEBUG: recommend_from_reaction called ===")
    print(f"search_all_families: {kwargs.get('search_all_families', False)}")
    
    result = original_recommend(*args, **kwargs)
    
    print(f"\nResult detection metadata:")
    det = result.get('detection', {})
    print(f"  auto_family: {det.get('auto_family')}")
    print(f"  rule_family: {det.get('rule_family')}")
    print(f"  detected_family: {det.get('detected_family')}")
    print(f"  analysis_module_used: {det.get('analysis_module_used')}")
    
    rc = det.get('reactant_classification')
    if rc:
        print(f"\n  Reactant classification found:")
        print(f"    reaction_type: {rc.get('reaction_type')}")
        print(f"    confidence: {rc.get('confidence')}")
    
    print(f"\nResult family: {result.get('family')}")
    print(f"Result detected_family: {result.get('detected_family')}")
    print("=== END DEBUG ===\n")
    
    return result

rec_module.recommend_from_reaction = debug_recommend

# Now run the test
from chemtools import chem

print("Testing cross-family recommendation with debug output:")
print("=" * 80)

result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=5,
    search_all_families=True
)

print("\nFinal result:")
print(f"  Family: {result.get('family')}")
print(f"  Detected Family: {result.get('detected_family')}")
if 'detection' in result:
    print(f"  Detection Auto Family: {result['detection'].get('auto_family')}")
