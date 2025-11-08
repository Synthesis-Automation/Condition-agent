from chemtools.recommend import UnifiedRecommender
from chemtools.rule.analyzer import FeatureAnalyzer

# Test with acid chloride (not carboxylic acid) + aniline
reaction = 'ClC(=O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1'

analyzer = FeatureAnalyzer()
features = analyzer.analyze_reaction(reaction)

print("="*80)
print("Test: Acid Chloride + Aniline (No Carboxylic Acid)")
print("="*80)
print(f"\nReaction: {reaction}")
print("\nDetected Features:")
print(f"  carboxylic_acid_present: {features.get('carboxylic_acid_present', False)}")
print(f"  aniline_present: {features.get('aniline_present', False)}")
print(f"  primary_amine_present: {features.get('primary_amine_present', False)}")

print("\nExpected: Amide formation rule should be FILTERED OUT")
print("Reason: applies_if requires 'carboxylic_acid_present' (ALL condition)")

rec = UnifiedRecommender()

# Without validation
results_no_val = rec.recommend(reaction, top_k=5, validate_rules=False)
amide_no = [r for r in results_no_val if r.id == 'amide_formation_v2']

print("\n" + "-"*80)
print("WITHOUT validation:")
if amide_no:
    print(f"  ✅ amide_formation_v2 found (similarity: {amide_no[0].similarity:.3f})")
else:
    print(f"  ℹ️  amide_formation_v2 not in top 5 (low DRFP similarity)")

# With validation
results_val = rec.recommend(reaction, top_k=5, validate_rules=True)
amide_yes = [r for r in results_val if r.id == 'amide_formation_v2']

print("\nWITH validation:")
if amide_yes:
    print(f"  ❌ amide_formation_v2 found (similarity: {amide_yes[0].similarity:.3f})")
    print(f"     ERROR: Should have been filtered out!")
else:
    print(f"  ✅ amide_formation_v2 correctly filtered out")

print("\n" + "="*80)
print("Top 5 results (with validation):")
print("="*80)
for r in results_val:
    print(f"{r.rank}. {r.name}")
    print(f"   Type: {r.source_type}, Family: {r.family}, Similarity: {r.similarity:.3f}")
