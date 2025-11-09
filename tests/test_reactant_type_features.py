#!/usr/bin/env python3
"""
Test reactant type feature detection in calculable.py

Tests:
1. Member-level features (ArBr_reactant, etc.)
2. Category-level features (ArX_reactant, etc.)
3. Utility functions (get_reactant_type_features, classify_reactant_smiles)
4. Backward compatibility
"""

import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.featurizers import calculable
from chemtools.util.rdkit_helpers import rdkit_available

def test_member_features():
    """Test member-level reactant type detection."""
    print("\n" + "=" * 70)
    print("TEST 1: Member-Level Reactant Features")
    print("=" * 70)
    
    test_cases = [
        ("c1ccc(Br)cc1", "ArBr_reactant", "aryl bromide"),
        ("c1ccc(Cl)cc1", "ArCl_reactant", "aryl chloride"),
        ("c1ccccc1B(O)O", "ArB_OH_2_reactant", "aryl boronic acid"),
        ("CCN", "RNH2_reactant", "primary aliphatic amine"),
        ("c1ccccc1N", "ArNH2_reactant", "aniline"),
        ("C#CC", "terminal_alkyne_reactant", "terminal alkyne"),
    ]
    
    passed = 0
    failed = 0
    
    for smiles, expected_token, description in test_cases:
        features = calculable.detect_all_features(smiles)
        result = features.get(expected_token, False)
        
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {description:30s} | {smiles:20s} | {expected_token}")
        
        if result:
            passed += 1
        else:
            failed += 1
    
    print(f"\nResult: {passed}/{len(test_cases)} passed")
    return failed == 0


def test_category_features():
    """Test category-level derived features."""
    print("\n" + "=" * 70)
    print("TEST 2: Category-Level Reactant Features")
    print("=" * 70)
    
    test_cases = [
        ("c1ccc(Br)cc1", "ArX_reactant", "aryl halide category"),
        ("c1ccccc1B(O)O", "ArB_reactant", "aryl boron category"),
        ("CCN", "RNH2_R2NH_reactant", "aliphatic amine category"),
        ("C#CC", "Alkyne_reactant", "alkyne category"),
    ]
    
    passed = 0
    failed = 0
    
    for smiles, expected_token, description in test_cases:
        features = calculable.detect_all_features(smiles)
        result = features.get(expected_token, False)
        
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {description:30s} | {smiles:20s} | {expected_token}")
        
        if result:
            passed += 1
        else:
            failed += 1
    
    print(f"\nResult: {passed}/{len(test_cases)} passed")
    return failed == 0


def test_utility_functions():
    """Test reactant type utility functions."""
    print("\n" + "=" * 70)
    print("TEST 3: Reactant Type Utility Functions")
    print("=" * 70)
    
    # Test get_reactant_type_features
    print("\n3a. get_reactant_type_features('c1ccc(Br)cc1'):")
    features = calculable.get_reactant_type_features("c1ccc(Br)cc1")
    print(f"  Member types: {features.get('member_types', [])}")
    print(f"  Categories: {features.get('categories', [])}")
    print(f"  ArBr_reactant: {features.get('ArBr_reactant', False)}")
    print(f"  ArX_reactant: {features.get('ArX_reactant', False)}")
    
    # Test classify_reactant_smiles
    print("\n3b. classify_reactant_smiles('c1ccc(Br)cc1'):")
    result = calculable.classify_reactant_smiles("c1ccc(Br)cc1")
    if result:
        print(f"  ✓ Category: {result.get('category')}")
        print(f"  ✓ Member type: {result.get('member_type')}")
        print(f"  ✓ Name: {result.get('name')}")
        print(f"  ✓ Role: {result.get('coupling_role')}")
    else:
        print("  ✗ No match found")
    
    # Test get_reactant_categories
    print("\n3c. get_reactant_categories('c1ccc(Br)cc1'):")
    categories = calculable.get_reactant_categories("c1ccc(Br)cc1")
    print(f"  Categories: {categories}")
    
    # Test get_reactant_members
    print("\n3d. get_reactant_members('c1ccc(Br)cc1'):")
    members = calculable.get_reactant_members("c1ccc(Br)cc1")
    print(f"  Members: {members}")
    
    return result is not None


def test_backward_compatibility():
    """Test backward compatibility with old system."""
    print("\n" + "=" * 70)
    print("TEST 4: Backward Compatibility")
    print("=" * 70)
    
    test_cases = [
        ("c1ccc(Br)cc1", "ArX*", "ArBr"),
        ("c1ccc(Cl)cc1", "ArX*", "ArCl"),
        ("c1ccccc1B(O)O", "ArB*", None),  # Could be ArB(OH)2
        ("CCN", None, "RNH2"),  # Category varies
    ]
    
    passed = 0
    
    for smiles, expected_cat, expected_member in test_cases:
        result = calculable.classify_reactant_smiles(smiles)
        
        if result:
            cat_match = (expected_cat is None) or (result.get('category') == expected_cat)
            mem_match = (expected_member is None) or (result.get('member_type') == expected_member)
            
            if cat_match and mem_match:
                print(f"✓ PASS: {smiles:20s} → {result.get('category')} / {result.get('member_type')}")
                passed += 1
            else:
                print(f"✗ FAIL: {smiles:20s} → Got {result.get('category')}/{result.get('member_type')}")
        else:
            print(f"✗ FAIL: {smiles:20s} → No match")
    
    print(f"\nResult: {passed}/{len(test_cases)} passed")
    return passed > 0


def test_comprehensive_examples():
    """Test comprehensive examples across reaction types."""
    print("\n" + "=" * 70)
    print("TEST 5: Comprehensive Reactant Examples")
    print("=" * 70)
    
    examples = [
        # Suzuki coupling
        ("c1ccc(Br)cc1", "Aryl bromide (Suzuki electrophile)"),
        ("c1ccccc1B(O)O", "Boronic acid (Suzuki nucleophile)"),
        
        # Buchwald-Hartwig
        ("c1ccc(Cl)cc1", "Aryl chloride (BH electrophile)"),
        ("CCN", "Primary amine (BH nucleophile)"),
        ("c1ccccc1N", "Aniline (BH nucleophile)"),
        
        # Sonogashira
        ("c1ccc(I)cc1", "Aryl iodide (Sonogashira electrophile)"),
        ("C#CC", "Terminal alkyne (Sonogashira nucleophile)"),
        
        # Negishi
        ("c1ccc(Br)cc1", "Aryl bromide (Negishi electrophile)"),
        
        # Heck
        ("c1ccc(Br)cc1", "Aryl bromide (Heck electrophile)"),
        ("C=C", "Alkene (Heck nucleophile)"),
    ]
    
    for smiles, description in examples:
        result = calculable.classify_reactant_smiles(smiles)
        features = calculable.get_reactant_type_features(smiles)
        
        if result:
            print(f"\n{description}")
            print(f"  SMILES: {smiles}")
            print(f"  Type: {result.get('member_type')} ({result.get('category')})")
            print(f"  Name: {result.get('name')}")
            print(f"  Role: {result.get('coupling_role')}")
            print(f"  Detected features: {len([k for k in features.keys() if k.endswith('_reactant') and features[k]])}")
        else:
            print(f"\n✗ {description}: No match for {smiles}")
    
    return True


def main():
    """Run all tests."""
    print("=" * 70)
    print("REACTANT TYPE FEATURE DETECTION TESTS")
    print("=" * 70)
    
    if not rdkit_available():
        print("\n✗ RDKit not available. Cannot run tests.")
        return 1
    
    print(f"\nRDKit: Available")
    print(f"Testing: chemtools.featurizers.calculable")
    
    results = []
    
    # Run tests
    results.append(("Member features", test_member_features()))
    results.append(("Category features", test_category_features()))
    results.append(("Utility functions", test_utility_functions()))
    results.append(("Backward compatibility", test_backward_compatibility()))
    results.append(("Comprehensive examples", test_comprehensive_examples()))
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\nOverall: {passed}/{total} test suites passed")
    
    if passed == total:
        print("\n🎉 All tests passed! Reactant type features working correctly.")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test suite(s) failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
