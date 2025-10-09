"""
Test filtering of precedents with unknown reagents (not in database).
"""

from chemtools.recommend.core import recommend_from_reaction
from chemtools.ml.simple_precedent_ranker import recommend_simple

# Test Suzuki reaction
suzuki_rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

print("=" * 80)
print("Testing filter_unknown_reagents parameter")
print("=" * 80)

print("\n1. recommend_from_reaction() WITHOUT filtering:")
print("-" * 80)
results_no_filter = recommend_from_reaction(
    reaction=suzuki_rxn,
    k=30,
    rerank_strategy='none',
    filter_unknown_reagents=False
)
precs_no_filter = results_no_filter.get('precedent_pack', {}).get('precedents', [])
print(f"Total precedents: {len(precs_no_filter)}")
print(f"Top 5 cores: {[p.get('core') for p in precs_no_filter[:5]]}")
print(f"Top 5 bases: {[p.get('base_uid') for p in precs_no_filter[:5]]}")
print(f"Top 5 solvents: {[p.get('solvent_uid') for p in precs_no_filter[:5]]}")

print("\n2. recommend_from_reaction() WITH filtering:")
print("-" * 80)
results_with_filter = recommend_from_reaction(
    reaction=suzuki_rxn,
    k=30,
    rerank_strategy='none',
    filter_unknown_reagents=True
)
precs_with_filter = results_with_filter.get('precedent_pack', {}).get('precedents', [])
print(f"Total precedents: {len(precs_with_filter)}")
print(f"Top 5 cores: {[p.get('core') for p in precs_with_filter[:5]]}")
print(f"Top 5 bases: {[p.get('base_uid') for p in precs_with_filter[:5]]}")
print(f"Top 5 solvents: {[p.get('solvent_uid') for p in precs_with_filter[:5]]}")
print(f"Filtered count: {len(precs_no_filter) - len(precs_with_filter)}")

print("\n3. recommend_simple() WITHOUT filtering:")
print("-" * 80)
simple_no_filter = recommend_simple(
    reaction_smiles=suzuki_rxn,
    family='Suzuki',
    k=30,
    rerank_strategy='none',
    filter_unknown_reagents=False
)
precs_simple_no = simple_no_filter.get('precedents', [])
print(f"Total precedents: {len(precs_simple_no)}")
print(f"Top cores: {simple_no_filter.get('top_cores', [])[:3]}")
print(f"Top bases: {simple_no_filter.get('top_bases', [])[:3]}")
print(f"Reasoning: {simple_no_filter.get('reasoning', [])}")

print("\n4. recommend_simple() WITH filtering:")
print("-" * 80)
simple_with_filter = recommend_simple(
    reaction_smiles=suzuki_rxn,
    family='Suzuki',
    k=30,
    rerank_strategy='none',
    filter_unknown_reagents=True
)
precs_simple_with = simple_with_filter.get('precedents', [])
print(f"Total precedents: {len(precs_simple_with)}")
print(f"Top cores: {simple_with_filter.get('top_cores', [])[:3]}")
print(f"Top bases: {simple_with_filter.get('top_bases', [])[:3]}")
print(f"Reasoning: {simple_with_filter.get('reasoning', [])}")
print(f"Filtered count: {len(precs_simple_no) - len(precs_simple_with)}")

print("\n" + "=" * 80)
print("✅ Test completed! Check if filtering reduced precedent count.")
print("=" * 80)
