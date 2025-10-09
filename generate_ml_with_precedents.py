"""
Generate a complete ML recommendation JSON file with precedents_used section
"""
import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from chemtools import chem

# Test Suzuki reaction
reaction = "Clc1ccc(C#N)cc1.c1coc(B(O)O)c1>>N#Cc1ccc(-c2ccoc2)cc1"

print("Generating ML recommendation with precedents_used...")
print(f"Reaction: {reaction}\n")

# Call the recommend.conditions function (same as CLI uses)
result = chem.recommend.conditions(
    reaction=reaction,
    reaction_type="Suzuki",
    k=25,
    limit=3,
    relax={},
    constraints={},
)

# Save to JSON file
output_file = Path("results/ml_recommendation_with_precedents_cli_test.json")
output_file.parent.mkdir(exist_ok=True)

with open(output_file, "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)

print(f"✅ JSON saved to: {output_file}")
print(f"\nFile size: {output_file.stat().st_size:,} bytes")

# Show structure
print("\nTop-level keys in output:")
for key in result.keys():
    print(f"  - {key}")

if "precedents_used" in result:
    prec_used = result["precedents_used"]
    print(f"\n✅ precedents_used section present!")
    print(f"   - total_count: {prec_used.get('total_count')}")
    print(f"   - top_precedents count: {len(prec_used.get('top_precedents', []))}")
    print(f"   - core_matched_precedents: {prec_used.get('core_matched_precedents', {}).get('count')} matching {prec_used.get('core_matched_precedents', {}).get('core')}")
    
    stats = prec_used.get("statistics", {})
    print(f"   - statistics: avg_yield={stats.get('average_yield')}%, range={stats.get('yield_range')}")
else:
    print("\n❌ precedents_used section NOT found!")

print("\nDone!")
