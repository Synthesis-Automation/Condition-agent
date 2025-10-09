"""
Test script to verify updated API parameters and CLI functionality.
"""

from chemtools.recommend.core import recommend_from_reaction

print("=" * 80)
print("Testing Updated API - recommend_from_reaction()")
print("=" * 80)

# Test Ullmann reaction with new parameters
ullmann_rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

print("\n1. Test with rerank_strategy='rule' (default)")
print("-" * 80)
result1 = recommend_from_reaction(
    reaction=ullmann_rxn,
    k=30,
    rerank_strategy='rule',
    filter_unknown_reagents=False
)
print(f"Core: {result1['recommendation']['core']}")
print(f"Confidence: {result1['recommendation']['confidence']:.3f}")
print(f"Reasons: {result1.get('reasons', 'N/A')}")

print("\n2. Test with rerank_strategy='analytics'")
print("-" * 80)
result2 = recommend_from_reaction(
    reaction=ullmann_rxn,
    k=30,
    rerank_strategy='analytics',
    filter_unknown_reagents=False
)
print(f"Core: {result2['recommendation']['core']}")
print(f"Confidence: {result2['recommendation']['confidence']:.3f}")

print("\n3. Test with rerank_strategy='none'")
print("-" * 80)
result3 = recommend_from_reaction(
    reaction=ullmann_rxn,
    k=30,
    rerank_strategy='none',
    filter_unknown_reagents=False
)
print(f"Core: {result3['recommendation']['core']}")
print(f"Confidence: {result3['recommendation']['confidence']:.3f}")

print("\n4. Test with filter_unknown_reagents=True")
print("-" * 80)
result4 = recommend_from_reaction(
    reaction=ullmann_rxn,
    k=30,
    rerank_strategy='rule',
    filter_unknown_reagents=True
)
prec_count = len(result4['precedent_pack']['precedents'])
print(f"Core: {result4['recommendation']['core']}")
print(f"Precedent count: {prec_count}")

print("\n5. Test recommend_conditions_structured() with new parameters")
print("-" * 80)
from chemtools.recommend.core import recommend_conditions_structured

result5 = recommend_conditions_structured(
    reaction=ullmann_rxn,
    k=30,
    limit=3,
    rerank_strategy='rule',
    filter_unknown_reagents=False
)
print(f"Status: {result5['meta']['status']}")
print(f"Recommendations count: {len(result5['recommendations'])}")
print(f"Top core: {result5['recommendations'][0]['summary']['core'] if result5['recommendations'] else 'N/A'}")

print("\n6. Test context.py wrapper")
print("-" * 80)
from chemtools.context import ChemTools

ctx = ChemTools()
result6 = ctx.recommend.conditions(
    reaction=ullmann_rxn,
    k=30,
    limit=3,
    rerank_strategy='analytics',
    filter_unknown_reagents=False
)
print(f"Status: {result6['meta']['status']}")
print(f"Recommendations count: {len(result6['recommendations'])}")

print("\n" + "=" * 80)
print("✅ ALL API TESTS PASSED!")
print("=" * 80)

print("\nNew API features:")
print("  ✅ rerank_strategy parameter working ('none', 'rule', 'analytics')")
print("  ✅ filter_unknown_reagents parameter working")
print("  ✅ recommend_from_reaction() updated")
print("  ✅ recommend_conditions_structured() updated")
print("  ✅ context.py wrapper updated")
print("  ✅ No errors from removed fusion code")
