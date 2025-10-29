"""
Test backwards compatibility of deprecated APIs.

Verifies that old functions still work and show deprecation warnings.
"""

import warnings
from chemtools.router import detect_family, detect_family_from_reaction
from chemtools.reaction_type_detector import detect_reaction_type


def test_detect_family_deprecated():
    """Test detect_family shows deprecation warning."""
    print("\n=== Test 1: detect_family() Deprecation ===")
    
    # Capture deprecation warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        result = detect_family(["Brc1ccccc1", "c1ccccc1B(O)O"])
        
        # Filter to only DeprecationWarnings (ignore taxonomy warnings)
        deprecation_warnings = [warning for warning in w if issubclass(warning.category, DeprecationWarning)]
        
        # Check warning was raised
        assert len(deprecation_warnings) == 1, f"Expected 1 DeprecationWarning, got {len(deprecation_warnings)}"
        assert "detect_family() is deprecated" in str(deprecation_warnings[0].message)
        assert "detect_reaction()" in str(deprecation_warnings[0].message)
        
        print(f"✓ Deprecation warning shown: {deprecation_warnings[0].message}")
    
    # Check function still works
    assert "family" in result
    assert "confidence" in result
    assert "hits" in result
    assert result["family"] == "suzuki_miyaura"
    
    print(f"✓ Function still works: family={result['family']}, conf={result['confidence']}")
    print("✓ PASSED")


def test_detect_family_from_reaction_deprecated():
    """Test detect_family_from_reaction shows deprecation warning."""
    print("\n=== Test 2: detect_family_from_reaction() Deprecation ===")
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        result = detect_family_from_reaction(
            "Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1",
            use_rxn_insight=False
        )
        
        # Filter to only DeprecationWarnings
        deprecation_warnings = [warning for warning in w if issubclass(warning.category, DeprecationWarning)]
        
        # Check warning
        assert len(deprecation_warnings) == 1, f"Expected 1 DeprecationWarning, got {len(deprecation_warnings)}"
        assert "detect_family_from_reaction() is deprecated" in str(deprecation_warnings[0].message)
        
        print(f"✓ Deprecation warning shown: {deprecation_warnings[0].message}")
    
    # Check function still works
    assert "family" in result
    assert "confidence" in result
    assert "hits" in result
    assert result["family"] == "suzuki_miyaura"
    
    # Check old schema fields present
    assert "rule" in result
    assert "status" in result
    
    print(f"✓ Function still works: family={result['family']}, conf={result['confidence']}")
    print(f"✓ Old schema preserved: rule={result['rule']['family']}, status={result['status']}")
    print("✓ PASSED")


def test_detect_reaction_type_deprecated():
    """Test detect_reaction_type shows deprecation warning."""
    print("\n=== Test 3: detect_reaction_type() Deprecation ===")
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        result = detect_reaction_type("Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1")
        
        # Filter to only DeprecationWarnings from our code (not from external libraries)
        deprecation_warnings = [
            warning for warning in w 
            if issubclass(warning.category, DeprecationWarning)
            and "detect_reaction_type() is deprecated" in str(warning.message)
        ]
        
        # Check warning
        assert len(deprecation_warnings) == 1, f"Expected 1 DeprecationWarning, got {len(deprecation_warnings)}"
        
        print(f"✓ Deprecation warning shown: {deprecation_warnings[0].message}")
    
    # Check function still works
    assert "available" in result
    assert "mapped_family" in result
    assert result["mapped_family"] == "suzuki_miyaura"
    
    # Check old schema fields
    assert "rxn_class" in result
    assert "rxn_name" in result
    assert "catalysts" in result
    
    print(f"✓ Function still works: family={result['mapped_family']}")
    print(f"✓ Old schema preserved: available={result['available']}")
    print("✓ PASSED")


def test_old_vs_new_api_equivalence():
    """Test that old and new APIs produce equivalent results."""
    print("\n=== Test 4: Old vs New API Equivalence ===")
    
    from chemtools import detect_reaction
    
    rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    # New API (no warnings)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        new_result = detect_reaction(rxn, use_ml=False)
        assert len(w) == 0, "New API should not raise warnings"
    
    # Old API (with warnings)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        old_result = detect_family_from_reaction(rxn, use_rxn_insight=False)
        assert len(w) == 1, "Old API should raise deprecation warning"
    
    # Compare results
    assert new_result["family"] == old_result["family"], \
        f"Family mismatch: new={new_result['family']}, old={old_result['family']}"
    
    assert abs(new_result["confidence"] - old_result["confidence"]) < 0.01, \
        f"Confidence mismatch: new={new_result['confidence']}, old={old_result['confidence']}"
    
    print(f"✓ Both APIs return: family={new_result['family']}, conf={new_result['confidence']:.2f}")
    print("✓ Results are equivalent")
    print("✓ PASSED")


def test_migration_path():
    """Test the migration path from old to new API."""
    print("\n=== Test 5: Migration Path ===")
    
    rxn = "Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1"
    
    # OLD WAY (deprecated)
    print("\nOLD (deprecated):")
    print("  from chemtools.router import detect_family_from_reaction")
    print(f"  result = detect_family_from_reaction('{rxn[:30]}...')")
    
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        old_result = detect_family_from_reaction(rxn, use_rxn_insight=False)
    
    print(f"  → family: {old_result['family']}")
    print(f"  → confidence: {old_result['confidence']}")
    
    # NEW WAY (recommended)
    print("\nNEW (recommended):")
    print("  from chemtools import detect_reaction")
    print(f"  result = detect_reaction('{rxn[:30]}...', use_ml=False)")
    
    from chemtools import detect_reaction
    new_result = detect_reaction(rxn, use_ml=False)
    
    print(f"  → family: {new_result['family']}")
    print(f"  → confidence: {new_result['confidence']}")
    
    print("\n✓ Migration path clear and documented")
    print("✓ PASSED")


if __name__ == "__main__":
    print("=" * 60)
    print("Backwards Compatibility & Deprecation Tests")
    print("=" * 60)
    
    try:
        test_detect_family_deprecated()
        test_detect_family_from_reaction_deprecated()
        test_detect_reaction_type_deprecated()
        test_old_vs_new_api_equivalence()
        test_migration_path()
        
        print("\n" + "=" * 60)
        print("✓ ALL BACKWARDS COMPATIBILITY TESTS PASSED!")
        print("=" * 60)
        print("\n📋 Summary:")
        print("  • Old APIs still work (backwards compatible)")
        print("  • Deprecation warnings shown correctly")
        print("  • Results match between old and new APIs")
        print("  • Migration path is clear")
        print("\n⚠️  Reminder: Old APIs will be removed in v2.0")
        
    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        raise
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        raise
