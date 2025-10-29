"""
Quick test script for Phase 1 detection implementation.

Tests the new detect_reaction() API with various reaction types.
"""

from chemtools import detect_reaction
import json


def test_suzuki():
    """Test Suzuki-Miyaura coupling detection."""
    print("\n=== Test 1: Suzuki-Miyaura Coupling ===")
    rxn = "Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1"
    result = detect_reaction(rxn, use_ml=False)
    
    print(f"Reaction: {rxn}")
    print(f"Family: {result['family']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Method: {result['method']}")
    print(f"Status: {result['status']}")
    
    assert result['family'] == 'suzuki_miyaura', f"Expected suzuki_miyaura, got {result['family']}"
    assert result['confidence'] >= 0.8, f"Confidence too low: {result['confidence']}"
    print("✓ PASSED")


def test_buchwald_hartwig():
    """Test Buchwald-Hartwig C-N coupling detection."""
    print("\n=== Test 2: Buchwald-Hartwig C-N Coupling ===")
    rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = detect_reaction(rxn, use_ml=False)
    
    print(f"Reaction: {rxn}")
    print(f"Family: {result['family']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Method: {result['method']}")
    print(f"Catalysts: {result['details']['catalysts']}")
    
    # Should detect cn_coupling (generic) without Pd catalyst in reaction
    # In real usage, catalyst would be in agents, boosting to buchwald_hartwig_c_n
    assert result['family'] in ['cn_coupling', 'buchwald_hartwig_c_n'], \
        f"Expected C-N coupling family, got {result['family']}"
    assert result['confidence'] >= 0.7, f"Confidence too low: {result['confidence']}"
    print("✓ PASSED")


def test_sonogashira():
    """Test Sonogashira coupling detection."""
    print("\n=== Test 3: Sonogashira Coupling ===")
    rxn = "Brc1ccccc1.C#C>>c1ccccc1C#C"
    result = detect_reaction(rxn, use_ml=False)
    
    print(f"Reaction: {rxn}")
    print(f"Family: {result['family']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Method: {result['method']}")
    
    assert result['family'] == 'sonogashira', f"Expected sonogashira, got {result['family']}"
    assert result['confidence'] >= 0.8, f"Confidence too low: {result['confidence']}"
    print("✓ PASSED")


def test_amide_coupling():
    """Test amide coupling detection."""
    print("\n=== Test 4: Amide Coupling ===")
    rxn = "CC(=O)O.CCN>>CC(=O)NCC"
    result = detect_reaction(rxn, use_ml=False)
    
    print(f"Reaction: {rxn}")
    print(f"Family: {result['family']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Method: {result['method']}")
    
    assert result['family'] == 'amide_coupling', f"Expected amide_coupling, got {result['family']}"
    assert result['confidence'] >= 0.7, f"Confidence too low: {result['confidence']}"
    print("✓ PASSED")


def test_grignard_addition():
    """Test Grignard addition to carbonyl."""
    print("\n=== Test 5: Grignard Addition ===")
    rxn = "CC(=O)C.CCMgBr>>CC(O)(CC)C"
    result = detect_reaction(rxn, use_ml=False)
    
    print(f"Reaction: {rxn}")
    print(f"Family: {result['family']}")
    print(f"Confidence: {result['confidence']}")
    print(f"Method: {result['method']}")
    
    assert result['family'] == 'grignard_addition', f"Expected grignard_addition, got {result['family']}"
    assert result['confidence'] >= 0.85, f"Confidence too low: {result['confidence']}"
    print("✓ PASSED")


def test_functional_groups():
    """Test functional group detection."""
    print("\n=== Test 6: Functional Group Detection ===")
    rxn = "Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1"
    result = detect_reaction(rxn, use_ml=False)
    
    fg = result['details']['functional_groups']
    print(f"Detected functional groups:")
    for group, present in fg.items():
        if present:
            print(f"  - {group}")
    
    assert fg['aryl_halide'], "Should detect aryl halide"
    assert fg['boron'], "Should detect boron"
    print("✓ PASSED")


def test_full_output():
    """Test complete output structure."""
    print("\n=== Test 7: Full Output Structure ===")
    rxn = "Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1"
    result = detect_reaction(rxn, use_ml=False)
    
    print("Full result structure:")
    print(json.dumps(result, indent=2, default=str))
    
    # Check required keys
    required_keys = ['family', 'confidence', 'method', 'status', 'details']
    for key in required_keys:
        assert key in result, f"Missing required key: {key}"
    
    detail_keys = ['reactants', 'catalysts', 'functional_groups', 'rule_prediction']
    for key in detail_keys:
        assert key in result['details'], f"Missing detail key: {key}"
    
    print("✓ PASSED")


if __name__ == "__main__":
    print("=" * 60)
    print("Phase 1 Detection System Tests")
    print("=" * 60)
    
    try:
        test_suzuki()
        test_buchwald_hartwig()
        test_sonogashira()
        test_amide_coupling()
        test_grignard_addition()
        test_functional_groups()
        test_full_output()
        
        print("\n" + "=" * 60)
        print("✓ ALL TESTS PASSED!")
        print("=" * 60)
        
    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        raise
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        raise
