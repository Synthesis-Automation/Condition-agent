#!/usr/bin/env python3
"""Debug the recommendation structure to see what's in the data."""

from chemtools import chem

result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=20,
    search_all_families=True
)

print("=" * 80)
print("INSPECTING RECOMMENDATION STRUCTURE")
print("=" * 80)

recommendations = result.get('recommendations', [])
print(f"\nTotal recommendations: {len(recommendations)}")

if recommendations:
    print("\n" + "=" * 80)
    print("FIRST RECOMMENDATION")
    print("=" * 80)
    
    rec = recommendations[0]
    
    print(f"\nKeys in recommendation: {list(rec.keys())}")
    
    print(f"\n--- summary ---")
    summary = rec.get('summary', {})
    for key, value in summary.items():
        print(f"  {key}: {value}")
    
    print(f"\n--- chemicals ---")
    chemicals = rec.get('chemicals', [])
    print(f"Total chemicals: {len(chemicals)}")
    for i, chem in enumerate(chemicals[:5], 1):  # First 5
        print(f"\n  Chemical {i}:")
        print(f"    name: {chem.get('name')}")
        print(f"    role: {chem.get('role')}")
        print(f"    abbreviation: {chem.get('abbreviation')}")
        print(f"    cas: {chem.get('cas')}")
    
    print(f"\n--- combo ---")
    combo = rec.get('combo', {})
    for key, value in combo.items():
        print(f"  {key}: {value}")
    
    print(f"\n--- conditions ---")
    conditions = rec.get('conditions', {})
    for key, value in conditions.items():
        print(f"  {key}: {value}")

# Check precedents
print("\n" + "=" * 80)
print("PRECEDENTS STRUCTURE")
print("=" * 80)

precedents_used = result.get('precedents_used', {})
top_precs = precedents_used.get('top_precedents', [])
print(f"\nTotal top precedents: {len(top_precs)}")

if top_precs:
    print(f"\n--- First Precedent ---")
    prec = top_precs[0]
    print(f"Keys: {list(prec.keys())}")
    print(f"\nreaction_id: {prec.get('reaction_id')}")
    print(f"catalyst_name: {prec.get('catalyst_name')}")
    print(f"base_name: {prec.get('base_name')}")
    print(f"solvent_name: {prec.get('solvent_name')}")
    print(f"condition_core: {prec.get('condition_core')}")
    
    # Check chemicals in precedent
    if 'chemicals' in prec:
        print(f"\nChemicals in precedent: {len(prec['chemicals'])}")
        for i, ch in enumerate(prec['chemicals'][:3], 1):
            print(f"  {i}. {ch.get('name')} - role: {ch.get('role')}")
