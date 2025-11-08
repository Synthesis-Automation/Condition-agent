from chemtools.recommend import UnifiedRecommender

rec = UnifiedRecommender()

# Without validation
results = rec.recommend(
    'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
    top_k=20,
    min_similarity=0.05,
    validate_rules=False
)

amide_rule = [r for r in results if r.id == 'amide_formation_v2']
print('Without validation:')
print(f'  amide_formation_v2 found: {len(amide_rule) > 0}')
if amide_rule:
    print(f'  Similarity: {amide_rule[0].similarity:.3f}')
    print(f'  Rank: {amide_rule[0].rank}')

print()

# With validation
results_v = rec.recommend(
    'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
    top_k=20,
    min_similarity=0.05,
    validate_rules=True
)

amide_rule_v = [r for r in results_v if r.id == 'amide_formation_v2']
print('With validation:')
print(f'  amide_formation_v2 found: {len(amide_rule_v) > 0}')
if amide_rule_v:
    print(f'  Similarity: {amide_rule_v[0].similarity:.3f}')
    print(f'  Rank: {amide_rule_v[0].rank}')
else:
    print('  ✅ Correctly filtered out!')
