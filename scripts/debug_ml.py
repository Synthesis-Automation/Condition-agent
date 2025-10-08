"""Quick test of ML recommendation system to debug empty results."""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from chemtools import recommend

# Test reaction
reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"

print("Testing ML Recommendation System")
print("="*60)
print(f"Reaction: {reaction}")
print("="*60)

# Call ML system
data = recommend.recommend_conditions_structured(
    reaction=reaction,
    reaction_type="C_N_Coupling_Pd",
    k=10,
    limit=6,
    relax={
        "use_drfp": False,
        "precompute_drfp": False,
        "use_rxn_insight": True,
    },
    constraints=None,
)

print("\nResult:")
print(f"Keys: {list(data.keys())}")
print(f"\nDetection: {data.get('detection', {})}")

# Check all possible keys
print(f"\nRecommendations count: {len(data.get('recommendations', []))}")
print(f"Alternatives count: {len(data.get('alternatives', []))}")
print(f"Precedents count: {len(data.get('precedents', []))}")

if data.get('recommendations'):
    print(f"\n✅ First recommendation:")
    import json
    print(json.dumps(data['recommendations'][0], indent=2, default=str))
else:
    print("\n⚠ NO RECOMMENDATIONS RETURNED!")
