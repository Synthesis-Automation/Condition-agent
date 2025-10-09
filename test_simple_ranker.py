"""
Test script to compare simple precedent ranking strategies.

Compares:
1. No reranking (similarity only)
2. Rule-based reranking
3. Analytics-based reranking
"""

from chemtools.ml.simple_precedent_ranker import recommend_simple

# Test reaction: Suzuki coupling
reaction_smiles = "Clc1ccc(C#N)cc1.c1coc(B(O)O)c1>>N#Cc1ccc(-c2ccoc2)cc1"
family = "Suzuki"

print("=" * 80)
print("SIMPLE PRECEDENT RANKING TEST")
print("=" * 80)
print(f"\nReaction: {reaction_smiles}")
print(f"Family: {family}\n")

# Strategy 1: No reranking (similarity only)
print("\n" + "=" * 80)
print("STRATEGY 1: SIMILARITY ONLY (no reranking)")
print("=" * 80)
result1 = recommend_simple(
    reaction_smiles=reaction_smiles,
    family=family,
    k=30,
    rerank_strategy='none'
)

print(f"\nMethod: {result1['method']}")
print("\nReasoning:")
for reason in result1['reasoning']:
    print(f"  - {reason}")

print("\nTop 5 Catalysts:")
for core, count in result1['top_cores']:
    print(f"  {count}x {core}")

print("\nTop 5 Bases:")
for base, count in result1['top_bases']:
    print(f"  {count}x CAS:{base}")

print("\nTop 5 Precedents (by similarity):")
for i, (prec, sim) in enumerate(zip(result1['precedents'][:5], result1['similarities'][:5]), 1):
    print(f"  {i}. Similarity: {sim:.3f} | Core: {prec.get('core')} | Base: {prec.get('base_uid')} | Yield: {prec.get('yield', 'N/A')}%")

# Strategy 2: Rule-based reranking
print("\n" + "=" * 80)
print("STRATEGY 2: RULE-BASED RERANKING")
print("=" * 80)
result2 = recommend_simple(
    reaction_smiles=reaction_smiles,
    family=family,
    k=30,
    rerank_strategy='rule'
)

print(f"\nMethod: {result2['method']}")
print("\nReasoning:")
for reason in result2['reasoning']:
    print(f"  - {reason}")

print("\nTop 5 Catalysts:")
for core, count in result2['top_cores']:
    print(f"  {count}x {core}")

print("\nTop 5 Bases:")
for base, count in result2['top_bases']:
    print(f"  {count}x CAS:{base}")

print("\nTop 5 Precedents (rule-reranked):")
for i, (prec, sim) in enumerate(zip(result2['precedents'][:5], result2['similarities'][:5]), 1):
    print(f"  {i}. Similarity: {sim:.3f} | Core: {prec.get('core')} | Base: {prec.get('base_uid')} | Yield: {prec.get('yield', 'N/A')}%")

# Strategy 3: Analytics-based reranking
print("\n" + "=" * 80)
print("STRATEGY 3: ANALYTICS-BASED RERANKING")
print("=" * 80)
result3 = recommend_simple(
    reaction_smiles=reaction_smiles,
    family=family,
    k=30,
    rerank_strategy='analytics'
)

print(f"\nMethod: {result3['method']}")
print("\nReasoning:")
for reason in result3['reasoning']:
    print(f"  - {reason}")

print("\nTop 5 Catalysts:")
for core, count in result3['top_cores']:
    print(f"  {count}x {core}")

print("\nTop 5 Bases:")
for base, count in result3['top_bases']:
    print(f"  {count}x CAS:{base}")

print("\nTop 5 Precedents (analytics-reranked):")
for i, (prec, sim) in enumerate(zip(result3['precedents'][:5], result3['similarities'][:5]), 1):
    print(f"  {i}. Similarity: {sim:.3f} | Core: {prec.get('core')} | Base: {prec.get('base_uid')} | Yield: {prec.get('yield', 'N/A')}%")

# Compare the three strategies
print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)

print("\n#1 Recommended Catalyst:")
print(f"  Similarity only:  {result1['top_cores'][0][0] if result1['top_cores'] else 'None'}")
print(f"  Rule-reranked:    {result2['top_cores'][0][0] if result2['top_cores'] else 'None'}")
print(f"  Analytics-reranked: {result3['top_cores'][0][0] if result3['top_cores'] else 'None'}")

print("\n#1 Recommended Base:")
print(f"  Similarity only:  {result1['top_bases'][0][0] if result1['top_bases'] else 'None'}")
print(f"  Rule-reranked:    {result2['top_bases'][0][0] if result2['top_bases'] else 'None'}")
print(f"  Analytics-reranked: {result3['top_bases'][0][0] if result3['top_bases'] else 'None'}")

print("\n" + "=" * 80)
print("RECOMMENDATION: Test both strategies on validation set to see which gives better results!")
print("=" * 80)
