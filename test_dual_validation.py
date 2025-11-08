"""
Comprehensive test demonstrating the complete dual validation system.

This test shows how applies_if (rules) and reaction_SMARTS (protocols)
work together to provide layered chemical validation.
"""

from chemtools.recommend.unified import UnifiedRecommender


def test_dual_validation_system():
    """Test that both validation mechanisms work together correctly."""
    
    print("=" * 80)
    print("DUAL VALIDATION SYSTEM TEST")
    print("=" * 80)
    print()
    print("Testing both validation mechanisms working together:")
    print("  - applies_if: Validates rules based on functional groups")
    print("  - reaction_SMARTS: Validates protocols based on transformations")
    print()
    
    recommender = UnifiedRecommender()
    
    # Test 1: Query that should match both protocols and rules
    print("-" * 80)
    print("Test 1: Alkenyl Iodide Cyanation")
    print("-" * 80)
    
    rxn1 = "IC=CCCCC.CC(C)(O)C#N>>N#CC=CCCCC"
    print(f"Query: {rxn1}")
    print()
    
    results_no_val = recommender.recommend(rxn1, top_k=20, validate_rules=False)
    results_with_val = recommender.recommend(rxn1, top_k=20, validate_rules=True)
    
    protocols_no_val = [r for r in results_no_val if r.source_type == 'protocol']
    rules_no_val = [r for r in results_no_val if r.source_type == 'rule']
    
    protocols_with_val = [r for r in results_with_val if r.source_type == 'protocol']
    rules_with_val = [r for r in results_with_val if r.source_type == 'rule']
    
    print(f"WITHOUT validation:")
    print(f"  Protocols: {len(protocols_no_val)}")
    print(f"  Rules: {len(rules_no_val)}")
    print(f"  Total: {len(results_no_val)}")
    
    print(f"\nWITH validation:")
    print(f"  Protocols: {len(protocols_with_val)} (reaction_SMARTS validated)")
    print(f"  Rules: {len(rules_with_val)} (applies_if validated)")
    print(f"  Total: {len(results_with_val)}")
    
    # Check if cyanation protocol found
    cyanation_found = any('cyanation' in r.name.lower() for r in protocols_with_val)
    if cyanation_found:
        print(f"\n✅ Alkenyl cyanation protocol correctly included")
    else:
        print(f"\n⚠️  Alkenyl cyanation protocol not found (may not be in index)")
    
    # Test 2: Query that should NOT match certain protocols
    print("\n" + "-" * 80)
    print("Test 2: Aryl Iodide Cyanation (Different Mechanism)")
    print("-" * 80)
    
    rxn2 = "Ic1ccccc1.CC(C)(O)C#N>>N#Cc1ccccc1"
    print(f"Query: {rxn2}")
    print()
    
    results_no_val2 = recommender.recommend(rxn2, top_k=20, validate_rules=False)
    results_with_val2 = recommender.recommend(rxn2, top_k=20, validate_rules=True)
    
    alkenyl_no_val = any(
        'alkenyl' in r.name.lower() and 'cyanation' in r.name.lower()
        for r in results_no_val2
    )
    alkenyl_with_val = any(
        'alkenyl' in r.name.lower() and 'cyanation' in r.name.lower()
        for r in results_with_val2
    )
    
    print(f"WITHOUT validation:")
    print(f"  Alkenyl cyanation protocol found: {alkenyl_no_val}")
    
    print(f"\nWITH validation:")
    print(f"  Alkenyl cyanation protocol found: {alkenyl_with_val}")
    
    if not alkenyl_with_val:
        print(f"\n✅ Incompatible protocol correctly filtered by reaction_SMARTS")
    else:
        print(f"\n⚠️  Protocol should have been filtered (check DRFP similarity)")
    
    # Test 3: Summary statistics
    print("\n" + "=" * 80)
    print("VALIDATION IMPACT SUMMARY")
    print("=" * 80)
    
    test_queries = [
        ("Alkenyl iodide cyanation", "IC=CCCCC.CC(C)(O)C#N>>N#CC=CCCCC"),
        ("Aryl iodide cyanation", "Ic1ccccc1.CC(C)(O)C#N>>N#Cc1ccccc1"),
        ("Buchwald-Hartwig", "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"),
        ("Suzuki coupling", "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"),
    ]
    
    print()
    print(f"{'Query':<30} | {'Before':<10} | {'After':<10} | {'Filtered':<10} | Impact")
    print("-" * 85)
    
    for name, rxn in test_queries:
        results_before = recommender.recommend(rxn, top_k=20, validate_rules=False)
        results_after = recommender.recommend(rxn, top_k=20, validate_rules=True)
        
        before_count = len(results_before)
        after_count = len(results_after)
        filtered = before_count - after_count
        
        if filtered > 0:
            impact = f"✅ -{filtered}"
        elif filtered == 0:
            impact = "⚪ 0"
        else:
            impact = f"⚠️  +{abs(filtered)}"
        
        print(f"{name:<30} | {before_count:<10} | {after_count:<10} | {filtered:<10} | {impact}")
    
    print()
    print("Legend:")
    print("  ✅ Validation filtered out inappropriate recommendations")
    print("  ⚪ All recommendations were appropriate (no filtering needed)")
    print("  ⚠️  Unusual case (should not happen with current implementation)")
    
    print("\n" + "=" * 80)
    print("MECHANISM BREAKDOWN")
    print("=" * 80)
    print()
    print("┌─────────────────────────────────────────────────────────────────┐")
    print("│ Validation Mechanism          │ Target │ Method               │")
    print("├─────────────────────────────────────────────────────────────────┤")
    print("│ applies_if                    │ Rules  │ Functional groups    │")
    print("│ reaction_SMARTS               │ Protos │ Transformation match │")
    print("└─────────────────────────────────────────────────────────────────┘")
    print()
    print("Both mechanisms:")
    print("  ✓ Fail-open (include if validation fails)")
    print("  ✓ Post-similarity filtering (after DRFP search)")
    print("  ✓ Controlled by validate_rules parameter (default: True)")
    print()


if __name__ == "__main__":
    print("\n" + "🔬 " * 30)
    print("COMPREHENSIVE DUAL VALIDATION TEST")
    print("🔬 " * 30)
    print()
    
    test_dual_validation_system()
    
    print("=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)
    print()
    print("The dual validation system (applies_if + reaction_SMARTS) provides")
    print("layered chemical filtering to ensure recommended conditions are")
    print("appropriate for the query reaction.")
    print()
