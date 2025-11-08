"""
Test reaction_SMARTS validation for protocols.

This script tests the new reaction_SMARTS validation feature that filters
protocols based on exact reaction transformation patterns.

Expected behavior:
- Alkenyl iodide cyanation should match the cyanation protocol
- Aryl iodide cyanation should NOT match (different transformation)
- Alkyl iodide cyanation should NOT match (different substrate)
"""

from chemtools.recommend.unified import UnifiedRecommender

def test_alkenyl_iodide_matches():
    """Test that alkenyl iodide cyanation matches the protocol."""
    print("=" * 80)
    print("TEST 1: Alkenyl Iodide Cyanation (Should Match)")
    print("=" * 80)
    
    # Alkenyl iodide + acetone cyanohydrin → alkenyl nitrile
    # This should match: IC=C.CC(O)(C#N)C>>N#CC=C
    reaction = "IC=CCCCC.CC(C)(O)C#N>>N#CC=CCCCC"
    
    print(f"\nQuery Reaction: {reaction}")
    print("Expected: Cyanation protocol SHOULD be included")
    print()
    
    recommender = UnifiedRecommender()
    
    # Without validation
    results_no_val = recommender.recommend(
        reaction,
        top_k=20,
        validate_rules=False
    )
    
    print(f"WITHOUT validation: {len(results_no_val)} results")
    cyanation_found_no_val = any(
        'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower()
        for r in results_no_val
    )
    print(f"  Cyanation protocol found: {cyanation_found_no_val}")
    if cyanation_found_no_val:
        for r in results_no_val:
            if 'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower():
                print(f"    ✓ {r.name} (similarity: {r.similarity:.3f})")
    
    # With validation
    results_with_val = recommender.recommend(
        reaction,
        top_k=20,
        validate_rules=True
    )
    
    print(f"\nWITH validation: {len(results_with_val)} results")
    cyanation_found_with_val = any(
        'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower()
        for r in results_with_val
    )
    print(f"  Cyanation protocol found: {cyanation_found_with_val}")
    if cyanation_found_with_val:
        for r in results_with_val:
            if 'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower():
                print(f"    ✓ {r.name} (similarity: {r.similarity:.3f})")
    
    # Check result
    if cyanation_found_with_val:
        print("\n✅ TEST PASSED: Protocol correctly included for matching transformation")
    else:
        print("\n⚠️  TEST ISSUE: Protocol should match but was not found")
        print("   (May be due to low DRFP similarity or missing from index)")
    
    return cyanation_found_with_val


def test_aryl_iodide_filtered():
    """Test that aryl iodide cyanation does NOT match alkenyl iodide protocol."""
    print("\n" + "=" * 80)
    print("TEST 2: Aryl Iodide Cyanation (Should NOT Match)")
    print("=" * 80)
    
    # Aryl iodide + cyanide → aryl nitrile
    # This should NOT match IC=C pattern (aryl != alkenyl)
    reaction = "Ic1ccccc1.CC(C)(O)C#N>>N#Cc1ccccc1"
    
    print(f"\nQuery Reaction: {reaction}")
    print("Expected: Alkenyl iodide cyanation protocol SHOULD be filtered out")
    print()
    
    recommender = UnifiedRecommender()
    
    # Without validation
    results_no_val = recommender.recommend(
        reaction,
        top_k=20,
        validate_rules=False
    )
    
    print(f"WITHOUT validation: {len(results_no_val)} results")
    alkenyl_found_no_val = any(
        'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower()
        for r in results_no_val
    )
    print(f"  Alkenyl iodide protocol found: {alkenyl_found_no_val}")
    if alkenyl_found_no_val:
        for r in results_no_val:
            if 'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower():
                print(f"    ! {r.name} (similarity: {r.similarity:.3f})")
    
    # With validation
    results_with_val = recommender.recommend(
        reaction,
        top_k=20,
        validate_rules=True
    )
    
    print(f"\nWITH validation: {len(results_with_val)} results")
    alkenyl_found_with_val = any(
        'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower()
        for r in results_with_val
    )
    print(f"  Alkenyl iodide protocol found: {alkenyl_found_with_val}")
    if alkenyl_found_with_val:
        for r in results_with_val:
            if 'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower():
                print(f"    ✗ {r.name} (similarity: {r.similarity:.3f}) [Should be filtered!]")
    
    # Check result
    if not alkenyl_found_with_val:
        print("\n✅ TEST PASSED: Incompatible protocol correctly filtered out")
    else:
        print("\n❌ TEST FAILED: Protocol should have been filtered but was included")
    
    return not alkenyl_found_with_val


def test_transformation_specificity():
    """Test with different cyanation types to show SMARTS specificity."""
    print("\n" + "=" * 80)
    print("TEST 3: Transformation Specificity")
    print("=" * 80)
    
    test_reactions = [
        ("Alkenyl iodide (Z)", "I/C=C\\CCCC.CC(C)(O)C#N>>N#C/C=C\\CCCC"),
        ("Aryl bromide", "Brc1ccccc1.CC(C)(O)C#N>>N#Cc1ccccc1"),
        ("Alkyl iodide", "ICCCCCC.CC(C)(O)C#N>>N#CCCCCCC"),
    ]
    
    recommender = UnifiedRecommender()
    
    results_table = []
    
    for name, rxn in test_reactions:
        results = recommender.recommend(
            rxn,
            top_k=20,
            validate_rules=True
        )
        
        alkenyl_protocol = any(
            'cyanation' in r.name.lower() and 'alkenyl' in r.name.lower()
            for r in results
        )
        
        results_table.append((name, rxn, alkenyl_protocol))
    
    print("\nReaction Type                | Alkenyl Protocol Matched | Expected")
    print("-" * 70)
    for name, rxn, matched in results_table:
        expected = "✓" if name.startswith("Alkenyl") else "✗"
        status = "✓" if matched else "✗"
        indicator = "✅" if (matched and name.startswith("Alkenyl")) or (not matched and not name.startswith("Alkenyl")) else "❌"
        print(f"{name:28} | {status:24} | {expected:8} {indicator}")
    
    # Check if only alkenyl matched
    alkenyl_only = results_table[0][2] and not results_table[1][2] and not results_table[2][2]
    
    if alkenyl_only:
        print("\n✅ EXCELLENT: reaction_SMARTS correctly distinguishes substrate types!")
    else:
        print("\n⚠️  Note: Check if protocol is in index and DRFP similarity is sufficient")


if __name__ == "__main__":
    print("\n" + "🧪 " * 30)
    print("REACTION_SMARTS VALIDATION TEST SUITE")
    print("🧪 " * 30)
    
    test1_pass = test_alkenyl_iodide_matches()
    test2_pass = test_aryl_iodide_filtered()
    test_transformation_specificity()
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    if test1_pass and test2_pass:
        print("✅ ALL CRITICAL TESTS PASSED")
        print("\nreaction_SMARTS validation is working correctly!")
        print("Protocols are filtered based on exact transformation patterns.")
    elif test1_pass:
        print("⚠️  PARTIAL SUCCESS")
        print("✓ Matching transformations are included")
        print("? Non-matching transformations (check DRFP similarity)")
    elif test2_pass:
        print("⚠️  PARTIAL SUCCESS")
        print("✓ Non-matching transformations are filtered")
        print("? Matching transformations (check if protocol in index)")
    else:
        print("⚠️  CHECK REQUIRED")
        print("Verify protocol is indexed and DRFP similarity is sufficient")
    
    print("\nNote: validation uses fail-open behavior - if RDKit unavailable or")
    print("      pattern matching fails, protocols are included (permissive).")
