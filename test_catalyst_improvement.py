"""
Test catalyst specificity improvement.

This test verifies that recommendations now show specific catalyst complexes
from precedents instead of generic "Palladium".
"""

from chemtools import recommend

# Suzuki reaction from the user's example
rxn = "C/C=C/Br.C=CB(O)O>>C/C=C/C=C"

print("Testing Catalyst Specificity Improvement")
print("=" * 80)
print(f"Reaction: {rxn}")
print()

# Run recommendation
result = recommend.recommend_from_reaction(
    reaction=rxn,
    k=50,
    family_override="Suzuki",
    max_variants=2,
    rerank_strategy='rule',
)

# Get formatted output
formatted = result.get('formatted', {})
recommendations = formatted.get('recommended_conditions', [])

print(f"Found {len(recommendations)} recommendations\n")

for rec in recommendations[:2]:
    rank = rec.get('rank', '?')
    chemicals = rec.get('chemicals', [])
    
    print(f"Rank {rank} Chemicals:")
    print("-" * 80)
    
    for chem in chemicals:
        role = chem.get('role', 'unknown')
        name = chem.get('name', 'N/A')
        cas = chem.get('cas', 'N/A')
        
        if role in ['metal_precursor', 'ligand', 'catalyst']:
            print(f"  ✓ {role.upper()}: {name}")
            print(f"    CAS: {cas}")
    print()

# Check precedents
precedents_used = formatted.get('precedents_used', {})
top_precedents = precedents_used.get('top_precedents', [])

if top_precedents:
    print("Top Precedent Catalyst Info:")
    print("-" * 80)
    top = top_precedents[0]
    catalysts = top.get('catalysts', [])
    for cat in catalysts:
        print(f"  Precedent catalyst: {cat.get('name')}")
        print(f"  CAS: {cat.get('cas')}")
        print(f"  Core: {top.get('core')}")
    print()

# Core matched precedents
core_matched = precedents_used.get('core_matched_precedents', {})
chosen_core = core_matched.get('core')
print(f"Chosen Core: {chosen_core}")
print(f"Core-matched precedents: {core_matched.get('count', 0)}")
print()

# Verify improvement
print("Verification:")
print("-" * 80)
if recommendations:
    rec1 = recommendations[0]
    catalysts = [c for c in rec1.get('chemicals', []) 
                 if c.get('role') in ['metal_precursor', 'ligand', 'catalyst']]
    
    if catalysts:
        for cat in catalysts:
            name = cat.get('name', '')
            if 'diphenylphosphino' in name.lower() or 'dppf' in name.lower():
                print("✅ SUCCESS: Specific catalyst complex shown in recommendation!")
                print(f"   Catalyst: {name}")
            elif name.lower() == 'palladium':
                print("⚠️  GENERIC: Still showing generic 'Palladium'")
                print(f"   This might be expected if no catalytic_system data in precedents")
            else:
                print(f"ℹ️  Catalyst: {name}")
    else:
        print("❌ No catalyst found in recommendation")
else:
    print("❌ No recommendations generated")
