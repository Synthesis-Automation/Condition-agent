#!/usr/bin/env python3
"""
Test Step 5: Rule-Based Reaction Family Detection

Tests:
- Rule-based reaction family detection (detect_family)
- SMARTS pattern matching for different reaction types
- Catalyst detection and override logic
- Confidence scoring for different reaction signatures
- Comparison between rule-based and ML-based detection (if available)

This validates the deterministic rule-based detection system as a fallback
when ML models are unavailable or for rapid classification.
"""

import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools import detect_reaction
from chemtools.recommend.core import recommend_from_reaction
from chemtools import reagent


def lookup_name(cas: str, role: str) -> str:
    """Look up chemical name from CAS number."""
    if not cas or cas == 'N/A':
        return 'N/A'
    
    info = reagent.enrich_reagent_info(cas, role)
    if info and info.get('name'):
        name = info['name']
        if name == cas or name.startswith('[Unknown'):
            return f"[Unknown {role}] {cas}"
        return f"{name} ({cas})"
    return f"[Unknown {role}] {cas}"


# Test reactions for different families
TEST_REACTIONS = {
    "Ullmann C-N (Cu)": {
        "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        "expected_family": "C_N_Coupling_Cu",  # After mapping
        "expected_rule_family": "Ullmann_CN",
        "confidence_min": 0.85,
        "description": "Aryl bromide + aniline �?Diphenylamine (Cu-catalyzed)"
    },
    "Buchwald-Hartwig C-N (Pd)": {
        "smiles": "Brc1ccccc1.Nc1ccccc1>[Pd]>c1ccccc1Nc1ccccc1",
        "expected_family": "C_N_Coupling_Pd",
        "expected_rule_family": "Buchwald_CN",
        "confidence_min": 0.85,
        "description": "Aryl bromide + aniline with Pd catalyst"
    },
    "Suzuki Coupling": {
        "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccccc1-c1ccccc1",
        "expected_family": "Suzuki_Miyaura",
        "expected_rule_family": "Suzuki_CC",
        "confidence_min": 0.85,
        "description": "Aryl bromide + boronic acid �?Biphenyl"
    },
    "Sonogashira Coupling": {
        "smiles": "Brc1ccccc1.C#C>>C#Cc1ccccc1",
        "expected_family": "Sonogashira",
        "expected_rule_family": "Sonogashira_CC",
        "confidence_min": 0.80,
        "description": "Aryl bromide + terminal alkyne �?Phenylacetylene"
    },
    "Amide Coupling": {
        "smiles": "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1",
        "expected_family": "Amide_Formation",
        "expected_rule_family": "Amide_Coupling",
        "confidence_min": 0.75,
        "description": "Acetic acid + benzylamine �?N-benzylacetamide"
    },
}


def test_5a_rule_based_detection():
    """Test 5a: Rule-based reaction family detection."""
    print("=" * 70)
    print("  Test Step 5: Rule-Based Reaction Family Detection")
    print("=" * 70)
    print()
    print("──" * 35)
    print("  Test 5a: Rule-Based Family Detection")
    print("──" * 35)
    print()
    
    passed_tests = 0
    total_tests = len(TEST_REACTIONS)
    
    for rxn_name, rxn_data in TEST_REACTIONS.items():
        print(f"   🧪 Testing: {rxn_name}")
        print(f"      Description: {rxn_data['description']}")
        print(f"      SMILES: {rxn_data['smiles'][:60]}{'...' if len(rxn_data['smiles']) > 60 else ''}")
        
        t_start = time.time()
        
        try:
            # Use detect_family_from_reaction for complete analysis
            result = detect_family_from_reaction(
                rxn_data['smiles'],
                use_rxn_insight=False  # Pure rule-based
            )
            
            t_elapsed = time.time() - t_start
            
            # Extract results
            detected_family = result.get('family', 'Unknown')
            confidence = result.get('confidence', 0.0)
            hits = result.get('hits', {})
            rule_info = result.get('rule', {})
            
            print(f"      ⏱️  Detection time: {t_elapsed:.3f}s")
            print(f"      📊 Results:")
            print(f"         Family: {detected_family}")
            print(f"         Confidence: {confidence:.2f}")
            
            # Show pattern hits
            pattern_hits = [k for k, v in hits.items() if v and not k.startswith('catalyst_')]
            if pattern_hits:
                print(f"         Pattern hits: {', '.join(pattern_hits)}")
            
            # Show catalyst detection
            catalysts = [k.replace('catalyst_', '').upper() for k, v in hits.items() if k.startswith('catalyst_') and v]
            if catalysts:
                print(f"         Catalysts detected: {', '.join(catalysts)}")
            
            # Validation
            expected_rule = rxn_data['expected_rule_family']
            min_confidence = rxn_data['confidence_min']
            
            family_match = detected_family == expected_rule
            confidence_ok = confidence >= min_confidence
            
            if family_match and confidence_ok:
                print(f"      �?PASSED")
                passed_tests += 1
            else:
                print(f"      �?FAILED")
                if not family_match:
                    print(f"         Expected family: {expected_rule}, got: {detected_family}")
                if not confidence_ok:
                    print(f"         Expected confidence >= {min_confidence}, got: {confidence:.2f}")
            
            print()
            
        except Exception as e:
            print(f"      �?ERROR: {e}")
            import traceback
            traceback.print_exc()
            print()
    
    print("──" * 35)
    print("  Summary")
    print("──" * 35)
    print()
    print(f"   Tests passed: {passed_tests}/{total_tests}")
    
    if passed_tests == total_tests:
        print("   �?Test 5a PASSED")
    else:
        print(f"   ⚠️  Test 5a PARTIAL ({passed_tests}/{total_tests} passed)")
    
    print()
    return 0 if passed_tests == total_tests else 1


def test_5b_pattern_recognition():
    """Test 5b: SMARTS pattern recognition."""
    print("──" * 35)
    print("  Test 5b: SMARTS Pattern Recognition")
    print("──" * 35)
    print()
    
    # Test specific patterns
    test_patterns = {
        "Aryl Halide": {
            "reactants": ["Brc1ccccc1"],
            "expected_hits": ["aryl_halide"],
            "not_expected": ["boron", "terminal_alkyne"]
        },
        "Boronic Acid": {
            "reactants": ["c1ccc(B(O)O)cc1"],
            "expected_hits": ["boron"],
            "not_expected": ["aryl_halide", "terminal_alkyne"]
        },
        "Terminal Alkyne": {
            "reactants": ["C#C"],
            "expected_hits": ["terminal_alkyne"],
            "not_expected": ["aryl_halide", "boron"]
        },
        "Aniline": {
            "reactants": ["Nc1ccccc1"],
            "expected_hits": ["nucleophile_n"],
            "not_expected": ["acid", "boron"]
        },
        "Carboxylic Acid": {
            "reactants": ["CC(=O)O"],
            "expected_hits": ["acid"],
            "not_expected": ["boron", "terminal_alkyne"]
        },
    }
    
    passed = 0
    total = len(test_patterns)
    
    for pattern_name, pattern_data in test_patterns.items():
        reactants = pattern_data['reactants']
        expected = set(pattern_data['expected_hits'])
        not_expected = set(pattern_data['not_expected'])
        
        print(f"   🔍 Testing pattern: {pattern_name}")
        print(f"      Reactants: {', '.join(reactants)}")
        
        result = detect_reaction(".".join(reactants) + ">>"", use_ml=False)
        hits = result.get('hits', {})
        
        # Check expected hits
        found_expected = {k for k in expected if hits.get(k)}
        found_unexpected = {k for k in not_expected if hits.get(k)}
        
        print(f"      Expected hits: {', '.join(expected)}")
        print(f"      Found: {', '.join(found_expected) if found_expected else 'none'}")
        
        if found_expected == expected and not found_unexpected:
            print(f"      �?PASSED")
            passed += 1
        else:
            print(f"      �?FAILED")
            if found_expected != expected:
                missing = expected - found_expected
                print(f"         Missing: {', '.join(missing)}")
            if found_unexpected:
                print(f"         Unexpected hits: {', '.join(found_unexpected)}")
        
        print()
    
    print("──" * 35)
    print("  Summary")
    print("──" * 35)
    print()
    print(f"   Pattern tests passed: {passed}/{total}")
    
    if passed == total:
        print("   �?Test 5b PASSED")
    else:
        print(f"   ⚠️  Test 5b PARTIAL ({passed}/{total} passed)")
    
    print()
    return 0 if passed == total else 1


def test_5c_catalyst_override():
    """Test 5c: Catalyst detection and override logic."""
    print("──" * 35)
    print("  Test 5c: Catalyst Detection & Override")
    print("──" * 35)
    print()
    
    # Test catalyst override for C-N coupling
    test_cases = [
        {
            "name": "Cu-catalyzed C-N (implicit)",
            "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "expected_family": "Ullmann_CN",
            "expected_catalyst": None,
            "description": "No explicit catalyst �?defaults to Ullmann"
        },
        {
            "name": "Cu-catalyzed C-N (explicit)",
            "smiles": "Brc1ccccc1.Nc1ccccc1>[Cu]>c1ccccc1Nc1ccccc1",
            "expected_family": "Ullmann_CN",
            "expected_catalyst": "Cu",
            "description": "Explicit Cu catalyst �?Ullmann_CN"
        },
        {
            "name": "Pd-catalyzed C-N (Buchwald-Hartwig)",
            "smiles": "Brc1ccccc1.Nc1ccccc1>[Pd]>c1ccccc1Nc1ccccc1",
            "expected_family": "Buchwald_CN",
            "expected_catalyst": "Pd",
            "description": "Explicit Pd catalyst �?Buchwald_CN override"
        },
    ]
    
    passed = 0
    total = len(test_cases)
    
    for test_case in test_cases:
        print(f"   🧪 {test_case['name']}")
        print(f"      {test_case['description']}")
        
        result = detect_family_from_reaction(
            test_case['smiles'],
            use_rxn_insight=False
        )
        
        detected_family = result.get('family', 'Unknown')
        catalysts = result.get('catalysts', {})
        hits = result.get('hits', {})
        
        print(f"      Detected family: {detected_family}")
        
        # Check catalyst detection from hits or catalysts dict
        detected_metals = []
        if catalysts and 'detected' in catalysts:
            detected_metals = catalysts.get('detected', [])
        else:
            # Check hits for catalyst_* flags
            detected_metals = [k.replace('catalyst_', '').upper() for k, v in hits.items() if k.startswith('catalyst_') and v]
        
        if detected_metals:
            print(f"      Detected catalysts: {', '.join(detected_metals)}")
        
        family_match = detected_family == test_case['expected_family']
        
        if test_case['expected_catalyst']:
            # Check both catalysts dict and hits
            catalyst_found = (
                test_case['expected_catalyst'] in detected_metals or
                test_case['expected_catalyst'].lower() in [m.lower() for m in detected_metals]
            )
        else:
            catalyst_found = True  # No catalyst expected
        
        if family_match and catalyst_found:
            print(f"      �?PASSED")
            passed += 1
        else:
            print(f"      �?FAILED")
            if not family_match:
                print(f"         Expected: {test_case['expected_family']}, Got: {detected_family}")
            if test_case['expected_catalyst'] and not catalyst_found:
                print(f"         Expected catalyst: {test_case['expected_catalyst']}")
        
        print()
    
    print("──" * 35)
    print("  Summary")
    print("──" * 35)
    print()
    print(f"   Catalyst tests passed: {passed}/{total}")
    
    if passed == total:
        print("   �?Test 5c PASSED")
    else:
        print(f"   ⚠️  Test 5c PARTIAL ({passed}/{total} passed)")
    
    print()
    return 0 if passed == total else 1


def test_5d_ml_vs_rule_comparison():
    """Test 5d: ML vs Rule-based detection comparison (if ML available)."""
    print("──" * 35)
    print("  Test 5d: ML vs Rule-Based Comparison")
    print("──" * 35)
    print()
    
    # Check if ML detection is available
    try:
        from chemtools._ml_helpers import is_available as rxn_avail
        ml_available = rxn_avail()
    except:
        ml_available = False
    
    if not ml_available:
        print("   ⚠️  ML detection (rxn-insight) not available")
        print("   ℹ️  Skipping ML comparison test")
        print("   �?Test 5d SKIPPED (not a failure)")
        print()
        return 0
    
    print("   �?ML detection available, comparing with rule-based...")
    print()
    
    # Test a few reactions with both methods
    test_rxn = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print(f"   🧪 Testing: {test_rxn}")
    print()
    
    # Rule-based only
    print("   📋 Rule-based detection:")
    t_start = time.time()
    result_rule = detect_reaction(test_rxn, use_ml=False)
    t_rule = time.time() - t_start
    
    print(f"      Family: {result_rule.get('family')}")
    print(f"      Confidence: {result_rule.get('confidence', 0):.2f}")
    print(f"      Time: {t_rule:.3f}s")
    print()
    
    # ML-based
    print("   🤖 ML-based detection:")
    t_start = time.time()
    result_ml = detect_reaction(test_rxn, use_ml=True)
    t_ml = time.time() - t_start
    
    print(f"      Family: {result_ml.get('family')}")
    print(f"      Confidence: {result_ml.get('confidence', 0):.2f}")
    print(f"      Time: {t_ml:.3f}s")
    
    if 'rxn' in result_ml:
        rxn_info = result_ml['rxn']
        print(f"      RXN Class: {rxn_info.get('rxn_class', 'N/A')}")
        print(f"      RXN Name: {rxn_info.get('rxn_name', 'N/A')}")
    
    print()
    
    # Compare
    agreement = result_rule.get('family') == result_ml.get('family')
    status = result_ml.get('status', 'unknown')
    
    print("   📊 Comparison:")
    print(f"      Agreement: {'�?Yes' if agreement else '⚠️  No'}")
    print(f"      Status: {status}")
    print(f"      Speed ratio: {t_ml/t_rule:.1f}x (ML is {t_ml/t_rule:.1f}x slower)")
    print()
    
    print("   �?Test 5d PASSED (informational)")
    print()
    return 0


def test_5e_detailed_conditions():
    """Test 5e: Get detailed reaction conditions with recommendations."""
    print("──" * 35)
    print("  Test 5e: Detailed Reaction Conditions")
    print("──" * 35)
    print()
    
    # Test a few reactions and get full recommendations
    test_cases = [
        {
            "name": "Ullmann C-N Coupling",
            "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        },
        {
            "name": "Suzuki Coupling",
            "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccccc1-c1ccccc1",
        },
    ]
    
    for test_case in test_cases:
        print(f"   🧪 {test_case['name']}")
        print(f"      SMILES: {test_case['smiles']}")
        print()
        
        # Step 1: Detect family
        print("      📋 Family Detection:")
        detection = detect_reaction(test_case['smiles'], use_ml=False)
        family = detection.get('family', 'Unknown')
        confidence = detection.get('confidence', 0.0)
        
        print(f"         Family: {family}")
        print(f"         Confidence: {confidence:.2f}")
        print()
        
        # Step 2: Get recommendations
        print("      💡 Recommended Conditions:")
        try:
            result = recommend_from_reaction(
                reaction=test_case['smiles'],
                k=10,
                max_variants=1,
            )
            
            # Extract recommendation
            if 'formatted' in result:
                formatted = result['formatted']
                conditions = formatted.get('recommended_conditions', [])
                
                if len(conditions) > 0:
                    cond = conditions[0]
                    combo = cond.get('combo', {})
                    chemicals = cond.get('chemicals', [])
                    temp_time = cond.get('conditions', {})
                    
                    # Show catalyst system
                    print("         🔬 Catalyst System:")
                    catalyst_items = [c for c in chemicals if c.get('role') in ['metal_precursor', 'ligand']]
                    if catalyst_items:
                        for item in catalyst_items:
                            role = item.get('role', 'unknown')
                            cas = item.get('cas', 'N/A')
                            name = lookup_name(cas, role) if cas != 'N/A' else 'N/A'
                            print(f"            {role.replace('_', ' ').title()}: {name}")
                    else:
                        core = combo.get('core', 'N/A')
                        print(f"            Core: {core}")
                    
                    # Show base
                    print()
                    print("         ⚗️  Base:")
                    base_items = [c for c in chemicals if c.get('role') == 'base']
                    if base_items:
                        for item in base_items:
                            cas = item.get('cas', 'N/A')
                            name = lookup_name(cas, 'base') if cas != 'N/A' else 'N/A'
                            print(f"            {name}")
                    else:
                        base_cas = combo.get('base_uid', 'N/A')
                        base_name = lookup_name(base_cas, 'base') if base_cas != 'N/A' else 'N/A'
                        print(f"            {base_name}")
                    
                    # Show solvent
                    print()
                    print("         🧪 Solvent:")
                    solvent_items = [c for c in chemicals if c.get('role') == 'solvent']
                    if solvent_items:
                        for item in solvent_items:
                            cas = item.get('cas', 'N/A')
                            name = lookup_name(cas, 'solvent') if cas != 'N/A' else 'N/A'
                            print(f"            {name}")
                    else:
                        solvent_cas = combo.get('solvent_uid', 'N/A')
                        solvent_name = lookup_name(solvent_cas, 'solvent') if solvent_cas != 'N/A' else 'N/A'
                        print(f"            {solvent_name}")
                    
                    # Show temperature and time
                    print()
                    print("         🌡�? Reaction Conditions:")
                    temperature = temp_time.get('temperature', 'N/A')
                    time = temp_time.get('time', 'N/A')
                    
                    # Also check combo for T_C and time_h
                    if temperature == 'N/A' or temperature is None:
                        temp_c = combo.get('T_C', 'N/A')
                        if temp_c and temp_c != 'N/A':
                            temperature = f"{temp_c}°C"
                    
                    if time == 'N/A' or time is None:
                        time_h = combo.get('time_h', 'N/A')
                        if time_h and time_h != 'N/A':
                            time = f"{time_h} h"
                    
                    print(f"            Temperature: {temperature}")
                    print(f"            Time: {time}")
                    
                    # Show precedent support
                    print()
                    print("         📊 Precedent Support:")
                    if 'precedent_pack' in result:
                        pack = result['precedent_pack']
                        support = pack.get('support', 0)
                        precedents = pack.get('precedents', [])
                        print(f"            Similar reactions: {support}")
                        print(f"            Top precedents: {len(precedents)}")
                        
                        # Show top precedent details
                        if precedents:
                            top = precedents[0]
                            score = top.get('score', 0)
                            print(f"            Best match score: {score:.3f}")
                else:
                    print("         ⚠️  No recommendations available")
            else:
                print("         ⚠️  No formatted output")
                
        except Exception as e:
            print(f"         �?Error: {e}")
        
        print()
        print("      " + "─" * 60)
        print()
    
    print("   �?Test 5e PASSED (informational)")
    print()
    return 0


def main():
    """Run all Test Step 5 tests."""
    print()
    
    # Test 5a: Rule-based detection
    result_5a = test_5a_rule_based_detection()
    
    # Test 5b: Pattern recognition
    result_5b = test_5b_pattern_recognition()
    
    # Test 5c: Catalyst override
    result_5c = test_5c_catalyst_override()
    
    # Test 5d: ML vs Rule comparison
    result_5d = test_5d_ml_vs_rule_comparison()
    
    # Test 5e: Detailed conditions
    result_5e = test_5e_detailed_conditions()
    
    # Overall summary
    print("=" * 70)
    print("  Test Step 5 Complete!")
    print("=" * 70)
    print()
    
    total_failed = result_5a + result_5b + result_5c
    
    if total_failed == 0:
        print("   🎉 All tests passed!")
    else:
        print(f"   ⚠️  Some tests failed or partially passed")
    
    print()
    return total_failed


if __name__ == "__main__":
    exit(main())

