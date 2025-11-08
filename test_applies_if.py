#!/usr/bin/env python
"""Test applies_if validation in UnifiedRecommender"""

from chemtools.recommend import UnifiedRecommender

def test_applies_if_validation():
    rec = UnifiedRecommender()
    
    print("="*80)
    print("Testing applies_if Validation")
    print("="*80)
    
    # Test 1: Valid amide formation (acid + amine)
    print("\n[Test 1] Valid Amide Formation Reaction")
    print("Query: O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1")
    print("Expected: Amide Formation rule should be included")
    
    results = rec.recommend(
        'O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1',
        top_k=5,
        validate_rules=True
    )
    
    amide_rules = [r for r in results if 'amide' in r.name.lower()]
    print(f"\nResult: Found {len(amide_rules)} amide rule(s)")
    if amide_rules:
        for r in amide_rules:
            print(f"  ✅ {r.name} (similarity: {r.similarity:.3f})")
    
    # Test 2: Invalid for amide formation (Suzuki - no acid or amine)
    print("\n" + "="*80)
    print("[Test 2] Invalid Amide Formation Reaction (Suzuki Coupling)")
    print("Query: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1")
    print("Expected: Amide Formation rule should be EXCLUDED")
    
    # Without validation
    results_no_val = rec.recommend(
        'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
        top_k=10,
        validate_rules=False
    )
    
    # With validation
    results_with_val = rec.recommend(
        'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
        top_k=10,
        validate_rules=True
    )
    
    amide_no_val = [r for r in results_no_val if 'amide' in r.name.lower()]
    amide_with_val = [r for r in results_with_val if 'amide' in r.name.lower()]
    
    print(f"\nWithout validation: {len(amide_no_val)} amide rule(s)")
    if amide_no_val:
        for r in amide_no_val:
            print(f"  ⚠️  {r.name} (similarity: {r.similarity:.3f})")
    
    print(f"\nWith validation: {len(amide_with_val)} amide rule(s)")
    if amide_with_val:
        for r in amide_with_val:
            print(f"  ❌ {r.name} (similarity: {r.similarity:.3f}) - Should be excluded!")
    else:
        print("  ✅ Correctly excluded amide formation rule")
    
    print("\nTop 3 validated results:")
    for r in results_with_val[:3]:
        print(f"  {r.rank}. {r.name} ({r.source_type}, sim: {r.similarity:.3f})")
    
    # Test 3: Check feature detection for Suzuki
    print("\n" + "="*80)
    print("[Test 3] Feature Detection for Suzuki Reaction")
    
    from chemtools.rule.analyzer import FeatureAnalyzer
    analyzer = FeatureAnalyzer()
    
    features = analyzer.analyze_reaction('Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1')
    
    print("\nDetected features:")
    print(f"  carboxylic_acid_present: {features.get('carboxylic_acid_present', False)}")
    print(f"  primary_amine_present: {features.get('primary_amine_present', False)}")
    print(f"  secondary_amine_present: {features.get('secondary_amine_present', False)}")
    print(f"  aniline_present: {features.get('aniline_present', False)}")
    print(f"  sp2_halide_present: {features.get('sp2_halide_present', False)}")
    print(f"  sp2_boron_present: {features.get('sp2_boron_present', False)}")
    
    print("\n" + "="*80)
    print("✅ Testing Complete")
    print("="*80)

if __name__ == '__main__':
    test_applies_if_validation()
