"""
Test simple ranker on Ullmann C-N coupling to verify Cu catalyst is preferred.

This tests that rule-based reranking correctly identifies Cu-based conditions
(not Pd-based) for Ullmann reactions.
"""

from chemtools.ml.simple_precedent_ranker import recommend_simple

# Ullmann C-N coupling: aryl halide + aniline
reaction_smiles = "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1"
family = "Ullmann_CN"

print("=" * 80)
print("ULLMANN C-N COUPLING TEST")
print("=" * 80)
print(f"\nReaction: {reaction_smiles}")
print(f"Family: {family}")
print("\nExpected: Cu-based catalyst (not Pd)")
print("=" * 80)

# Test all three strategies
strategies = ['none', 'rule', 'analytics']

results = {}
for strategy in strategies:
    print(f"\n{'='*80}")
    print(f"STRATEGY: {strategy.upper()}")
    print('='*80)
    
    result = recommend_simple(
        reaction_smiles=reaction_smiles,
        family=family,
        k=30,
        rerank_strategy=strategy
    )
    
    results[strategy] = result
    
    print(f"\nReasoning:")
    for reason in result['reasoning']:
        print(f"  - {reason}")
    
    print(f"\nTop 5 Catalysts:")
    for core, count in result['top_cores']:
        # Check if Cu or Pd
        is_cu = 'Cu' in core or 'cu' in core.lower()
        is_pd = 'Pd' in core or 'pd' in core.lower()
        marker = "✅ CORRECT (Cu)" if is_cu else ("❌ WRONG (Pd)" if is_pd else "")
        print(f"  {count}x {core} {marker}")
    
    print(f"\nTop 3 Precedents:")
    for i, (prec, sim) in enumerate(zip(result['precedents'][:3], result['similarities'][:3]), 1):
        core = prec.get('core', 'Unknown')
        is_cu = 'Cu' in core or 'cu' in core.lower()
        is_pd = 'Pd' in core or 'pd' in core.lower()
        marker = "✅" if is_cu else ("❌" if is_pd else "❓")
        print(f"  {i}. {marker} Sim: {sim:.3f} | Core: {core} | Yield: {prec.get('yield', 'N/A')}%")

# Final comparison
print("\n" + "=" * 80)
print("FINAL COMPARISON")
print("=" * 80)

for strategy in strategies:
    result = results[strategy]
    if result['top_cores']:
        top_core = result['top_cores'][0][0]
        is_cu = 'Cu' in top_core or 'cu' in top_core.lower()
        status = "✅ CORRECT" if is_cu else "❌ WRONG"
        print(f"\n{strategy.upper():15s} → {top_core:40s} {status}")
    else:
        print(f"\n{strategy.upper():15s} → No catalysts found")

print("\n" + "=" * 80)
print("EXPECTED: Rule-based reranking should boost Cu catalysts to top")
print("=" * 80)
