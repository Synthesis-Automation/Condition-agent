"""
Test that precedents_used is included in ML recommendation output via the CLI path
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools import chem

# Test Suzuki reaction
reaction = "Clc1ccc(C#N)cc1.c1coc(B(O)O)c1>>N#Cc1ccc(-c2ccoc2)cc1"

print("Testing ML recommendation with precedents_used...")
print(f"Reaction: {reaction}")
print()

# Call the recommend.conditions function (same as CLI uses)
result = chem.recommend.conditions(
    reaction=reaction,
    reaction_type="Suzuki",
    k=25,
    limit=3,
    relax={},
    constraints={},
)

# Check if precedents_used is in the output
if "precedents_used" in result:
    precedents_used = result["precedents_used"]
    print("✅ precedents_used section found in output!")
    print(f"   Total precedents: {precedents_used.get('total_count', 'N/A')}")
    
    top_precs = precedents_used.get("top_precedents", [])
    print(f"   Top precedents count: {len(top_precs)}")
    
    if top_precs:
        print("\n   First precedent:")
        first = top_precs[0]
        print(f"   - Reaction ID: {first.get('reaction_id')}")
        print(f"   - Core: {first.get('core')}")
        print(f"   - Yield: {first.get('yield')}%")
        
        if first.get('catalysts'):
            print(f"   - Catalysts: {[c.get('name') for c in first.get('catalysts', [])]}")
        if first.get('reagents'):
            print(f"   - Reagents: {[r.get('name') for r in first.get('reagents', [])]}")
        if first.get('solvents'):
            print(f"   - Solvents: {[s.get('name') for s in first.get('solvents', [])]}")
    
    stats = precedents_used.get("statistics", {})
    if stats:
        print(f"\n   Statistics:")
        print(f"   - Average yield: {stats.get('average_yield', 'N/A')}%")
        print(f"   - Yield range: {stats.get('yield_range', 'N/A')}")
    
    print("\n✅ Test PASSED - precedents_used is properly included!")
else:
    print("❌ precedents_used section NOT found in output")
    print("\nAvailable keys:", list(result.keys()))
    print("\n❌ Test FAILED")
