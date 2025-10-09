"""
Test to verify Ullmann reactions get Cu-based conditions, not Suzuki/Pd.
"""

import json
from chemtools.ml.fusion_recommender import recommend_with_fusion

# Test Ullmann C-N coupling: bromopyrimidine + aniline
ullmann_reaction = "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1"

print("Testing Ullmann C-N coupling reaction...")
print(f"Reaction: {ullmann_reaction}\n")

# Run fusion recommendation
result = recommend_with_fusion(
    reaction_smiles=ullmann_reaction,
    family="Ullmann_CN",
    k=50,
    top_n=3
)

print("=" * 60)
print("FUSION RESULTS")
print("=" * 60)

# Check detection - it's in the evidence
evidence = result.get('evidence', {})
detection = evidence.get('detection', {}) if isinstance(evidence, dict) else {}
print(f"\nDetected family: {detection.get('family', 'N/A')}")
print(f"Auto-detected: {detection.get('auto_family', 'N/A')}")
print(f"Rule family: {detection.get('rule_family', 'N/A')}")

# Check recommendations - they're at the top level!
recommendations = result.get('recommended_conditions', [])

print(f"\nResult keys: {list(result.keys())}")
print(f"Number of recommendations: {len(recommendations)}")

if not recommendations:
    print("\n[WARNING] No recommendations found in result!")
    print("\nFull result structure:")
    print(json.dumps(result, indent=2)[:2000])
    print("...(truncated)")

if recommendations:
    top_rec = recommendations[0]
    print(f"\n{'=' * 60}")
    print("TOP RECOMMENDATION (Rank 1)")
    print(f"{'=' * 60}")
    
    # Extract from ScoredCandidate.candidate
    cand = top_rec.candidate
    print(f"Core: {cand.core}")
    print(f"Base: {cand.base}")
    print(f"Solvent: {cand.solvent}")
    
    # Check if it's Cu-based (correct) or Pd-based (wrong)
    core_name = str(cand.core or '')
    if 'Cu' in core_name or 'Copper' in core_name:
        print("\n[SUCCESS] Cu-based catalyst - Ullmann reaction (CORRECT!)")
    elif 'Pd' in core_name or 'Palladium' in core_name:
        print("\n[FAIL] Pd-based catalyst - Suzuki reaction (BUG NOT FIXED!)")
    else:
        print(f"\n[UNKNOWN] Catalyst type: {core_name}")
    
    if hasattr(cand, 'T_C') and cand.T_C is not None:
        print(f"Temperature: {cand.T_C} C")
    
    if hasattr(cand, 'time_h') and cand.time_h is not None:
        print(f"Time: {cand.time_h} h")
    
    print(f"\nFusion score: {top_rec.total_score:.3f}")
    print(f"Confidence: {top_rec.confidence}")

# Check fusion meta - it's at top level
adaptive_weights = result.get('adaptive_weights', {})
reasoning = result.get('reasoning', [])

if reasoning:
    print(f"\n{'=' * 60}")
    print("FUSION REASONING")
    print(f"{'=' * 60}")
    for reason in reasoning:
        print(f"  • {reason}")

weights = adaptive_weights
if weights:
    print(f"\n{'=' * 60}")
    print("ADAPTIVE WEIGHTS")
    print(f"{'=' * 60}")
    print(f"  Precedents (α): {weights.get('α', 0):.3f}")
    print(f"  Analytics (β):  {weights.get('β', 0):.3f}")
    print(f"  Rules (γ):      {weights.get('γ', 0):.3f}")
    print(f"  ML (δ):         {weights.get('δ', 0):.3f}")

print("\n" + "=" * 60)
print("TEST COMPLETE")
print("=" * 60)
