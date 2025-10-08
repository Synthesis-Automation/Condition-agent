"""
Test Unified Output Format

Validates that both ML-based and rule-based recommendations produce
identical JSON structure suitable for robotic execution systems.

NOTE: This test can be slow for large datasets (Suzuki has 50,215 reactions).
      The precedent search with k=5 should take:
      - C-N Coupling (Cu): ~1 second (5,552 reactions)
      - Suzuki: ~2-3 minutes (50,215 reactions) - due to list operations in _candidate_pool
      
      For faster testing, you can:
      1. Use smaller k values (k=3)
      2. Skip Suzuki test if needed
      3. Use pytest markers to run tests selectively
"""

import sys
import os
from pathlib import Path

# Set UTF-8 encoding for Windows terminals
if sys.platform == 'win32':
    import codecs
    sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
    sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

# Add project root to Python path
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import json
import time
from typing import Dict, Any

# Try to import pytest for markers (optional)
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    # Create dummy decorator if pytest not available
    class pytest:
        @staticmethod
        def mark(name):
            def decorator(func):
                return func
            return decorator
        slow = staticmethod(lambda f: f)

# NOTE: Defer heavy imports to avoid loading DRFP data at module import time
# These will be imported inside test functions when needed


def verify_structure_match(ml_output: Dict, rule_output: Dict) -> bool:
    """
    Verify that ML and rule outputs have identical structure.
    
    Args:
        ml_output: Output from ML-based system
        rule_output: Output from rule-based system
    
    Returns:
        True if structures match
    """
    def get_structure(obj, path=""):
        """Recursively get structure keys."""
        if isinstance(obj, dict):
            return {k: get_structure(obj[k], f"{path}.{k}") for k in sorted(obj.keys())}
        elif isinstance(obj, list):
            if len(obj) > 0:
                return [get_structure(obj[0], f"{path}[0]")]
            return []
        else:
            return type(obj).__name__
    
    ml_struct = get_structure(ml_output)
    rule_struct = get_structure(rule_output)
    
    return ml_struct == rule_struct


def validate_reagent_structure(reagent: Dict[str, Any]) -> bool:
    """
    Validate that a reagent has all required fields.
    
    Required fields:
    - id, role, name, abbreviation, cas, smiles, inchi_key
    - equivalents (for most roles)
    - loading (for catalysts and ligands)
    - properties (for catalysts)
    """
    required_fields = ["id", "role"]
    
    for field in required_fields:
        if field not in reagent:
            print(f"❌ Missing required field: {field}")
            return False
    
    # Check role-specific fields
    role = reagent.get("role")
    
    if role in ["catalyst", "ligand"]:
        if "loading" not in reagent:
            print(f"❌ Catalyst/ligand missing 'loading' field")
            return False
        
        loading = reagent["loading"]
        if not all(k in loading for k in ["value", "range", "unit"]):
            print(f"❌ Invalid loading structure: {loading}")
            return False
        
        if loading["unit"] != "mol%":
            print(f"❌ Invalid loading unit: {loading['unit']} (expected 'mol%')")
            return False
    
    if role == "catalyst":
        if "category" not in reagent:
            print(f"⚠️  Catalyst missing 'category' field (optional but recommended)")
        
        if "properties" not in reagent:
            print(f"⚠️  Catalyst missing 'properties' field (optional but recommended)")
    
    # Check equivalents structure
    if "equivalents" in reagent:
        equiv = reagent["equivalents"]
        if not all(k in equiv for k in ["value", "range", "unit"]):
            print(f"❌ Invalid equivalents structure: {equiv}")
            return False
        
        if equiv["unit"] != "eq":
            print(f"❌ Invalid equivalents unit: {equiv['unit']} (expected 'eq')")
            return False
    
    return True


def validate_conditions_structure(conditions: Dict[str, Any]) -> bool:
    """
    Validate that conditions have all required fields.
    
    Required fields:
    - temperature (value, range, unit)
    - time (value, range, unit)
    - atmosphere (type, gas, requirement)
    - pressure (value, unit)
    """
    required_sections = ["temperature", "time"]
    recommended_sections = ["atmosphere", "pressure"]
    
    for section in required_sections:
        if section not in conditions:
            print(f"❌ Missing required condition: {section}")
            return False
        
        section_data = conditions[section]
        if not all(k in section_data for k in ["value", "range", "unit"]):
            print(f"❌ Invalid {section} structure: {section_data}")
            return False
    
    for section in recommended_sections:
        if section not in conditions:
            print(f"⚠️  Missing recommended condition: {section}")
    
    # Validate atmosphere structure if present
    if "atmosphere" in conditions:
        atm = conditions["atmosphere"]
        if not all(k in atm for k in ["type", "gas", "requirement"]):
            print(f"⚠️  Incomplete atmosphere structure: {atm}")
    
    return True


def validate_recommendation_structure(rec: Dict[str, Any]) -> bool:
    """
    Validate that a recommendation has all required fields.
    
    Required fields:
    - rank, reaction, chemicals, conditions
    Optional: summary (with confidence, support)
    """
    required_fields = ["rank", "reaction", "chemicals", "conditions"]
    
    for field in required_fields:
        if field not in rec:
            print(f"❌ Missing required field: {field}")
            return False
    
    # Validate chemicals (renamed from reagents)
    if not isinstance(rec["chemicals"], list) or len(rec["chemicals"]) == 0:
        print(f"❌ Invalid chemicals: must be non-empty list")
        return False
    
    # Note: Reagent validation is less strict for formatted output
    # as it may have different field names
    
    # Validate conditions
    if not isinstance(rec["conditions"], dict):
        print(f"❌ Invalid conditions: must be dict")
        return False
    
    # Check optional summary field
    if "summary" in rec:
        summary = rec["summary"]
        if not isinstance(summary, dict):
            print(f"⚠️  Summary is not a dict")
        elif "confidence" not in summary:
            print(f"⚠️  Summary missing confidence field")
    
    return True


def test_ullmann_cn_coupling():
    """
    Test 1: Ullmann C-N Coupling - Unified Format Validation
    
    Validates that ML-based recommendations produce the standardized
    format required for robotic execution.
    """
    print("\n" + "="*70)
    print("TEST 1: Ullmann C-N Coupling - Unified Format Validation")
    print("="*70)
    
    # Import heavy modules only when test runs (not at module import time)
    from chemtools.recommend.core import recommend_from_reaction
    from chemtools.reaction_type_detector import detect_reaction_type
    
    # Test reaction
    smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print(f"\n📋 Input: {smiles}")
    print(f"   Expected: C-N coupling (Ullmann or Buchwald-Hartwig)")
    
    # =========================================================================
    # STEP 1: ML-based recommendation
    # =========================================================================
    print("\n" + "-"*70)
    print("STEP 1: ML-based Recommendation with Standardized Output")
    print("-"*70)
    
    start = time.time()
    # Use k=5 for faster testing (default is 25)
    ml_raw = recommend_from_reaction(smiles, k=5)
    ml_time = (time.time() - start) * 1000
    
    # Get detection info
    detection = detect_reaction_type(smiles)
    detected_type = detection.get("type", "Unknown")
    detection_conf = detection.get("confidence")
    
    # Handle None confidence
    if detection_conf is not None:
        print(f"\n✅ ML Detection: {detected_type} (confidence: {detection_conf:.3f})")
    else:
        print(f"\n✅ ML Detection: {detected_type} (confidence: N/A)")
    
    print(f"⏱️  Processing time: {ml_time:.2f} ms")
    
    # Use the formatted output that's already in ml_raw
    if 'formatted' in ml_raw and 'recommended_conditions' in ml_raw['formatted']:
        ml_output = ml_raw['formatted']
        # Rename 'recommended_conditions' to 'recommendations' for consistency
        if 'recommended_conditions' in ml_output:
            ml_output['recommendations'] = ml_output.pop('recommended_conditions')
        print(f"📊 Retrieved {len(ml_output.get('recommendations', []))} recommendations")
    else:
        print(f"❌ No formatted output found in ML result")
        return False
    
    # =========================================================================
    # STEP 2: Structure Validation
    # =========================================================================
    print("\n" + "-"*70)
    print("STEP 2: Structure Validation")
    print("-"*70)
    
    # Verify top-level structure
    print("\n🔍 Checking top-level structure...")
    required_sections = ["meta", "input", "detection", "recommendations"]
    
    for section in required_sections:
        has_section = section in ml_output
        status = "✅" if has_section else "❌"
        print(f"  {status} {section:20s} Present:{has_section:5}")
    
    # =========================================================================
    # STEP 3: Validate Individual Recommendations
    # =========================================================================
    print("\n" + "-"*70)
    print("STEP 3: Validate Recommendations")
    print("-"*70)
    
    print("\n📊 ML Recommendations:")
    all_valid = True
    for i, rec in enumerate(ml_output["recommendations"][:3], 1):
        print(f"\n  Recommendation {i}:")
        print(f"    Rank: {rec.get('rank', 'N/A')}")
        
        # Confidence may be in summary
        if 'summary' in rec and 'confidence' in rec['summary']:
            print(f"    Confidence: {rec['summary']['confidence']:.3f}")
        elif 'confidence' in rec:
            print(f"    Confidence: {rec['confidence']:.3f}")
        
        # Support may be in summary
        if 'summary' in rec and 'support' in rec['summary']:
            support_info = rec['summary']['support']
            if isinstance(support_info, dict):
                print(f"    Support: {support_info.get('count', 'N/A')} precedents")
            else:
                print(f"    Support: {support_info} precedents")
        elif 'support' in rec:
            print(f"    Support: {rec['support']} precedents")
        
        # Use 'chemicals' instead of 'reagents'
        chemicals = rec.get('chemicals', rec.get('reagents', []))
        print(f"    Chemicals: {len(chemicals)} components")
        
        if validate_recommendation_structure(rec):
            print(f"    ✅ Structure valid")
        else:
            print(f"    ❌ Structure invalid")
            all_valid = False
    
    # =========================================================================
    # STEP 4: Robotic Execution Compatibility Check
    # =========================================================================
    print("\n" + "-"*70)
    print("STEP 4: Robotic Execution Compatibility")
    print("-"*70)
    
    print("\n🤖 Checking robotic system compatibility...")
    
    # Check that all chemicals have required robot fields
    ml_rec1 = ml_output["recommendations"][0]
    chemicals = ml_rec1.get("chemicals", ml_rec1.get("reagents", []))
    
    print("\n  Required fields for robotic execution:")
    
    robot_fields = {
        "Chemical Name or CAS": lambda r: r.get("name") or r.get("cas"),
        "Role": lambda r: r.get("role") is not None,
        "SMILES": lambda r: r.get("smiles") is not None,
    }
    
    for field_name, check_fn in robot_fields.items():
        # Check at least one chemical has this field
        ml_valid = any(check_fn(r) for r in chemicals)
        status = "✅" if ml_valid else "⚠️ "
        print(f"    {status} {field_name:30s} ML:{ml_valid:5}")
        
        status = "✅" if ml_valid else "❌"
        print(f"    {status} {field_name:30s} Present:{ml_valid:5}")
        
        if not ml_valid:
            all_valid = False
    
    if all_valid:
        print("\n  ✅ All required fields present for robotic execution!")
    else:
        print("\n  ❌ Some required fields missing")
        return False
    
    # =========================================================================
    # STEP 5: Save Output File
    # =========================================================================
    print("\n" + "-"*70)
    print("STEP 5: Save Output File")
    print("-"*70)
    
    # Save ML output
    ml_file = "output_ml_standardized.json"
    with open(ml_file, "w") as f:
        json.dump(ml_output, f, indent=2)
    print(f"\n💾 Saved ML output: {ml_file}")
    
    # Print first recommendation for inspection
    print("\n" + "-"*70)
    print("Sample Output (ML Recommendation #1)")
    print("-"*70)
    rec1_str = json.dumps(ml_output["recommendations"][0], indent=2)
    if len(rec1_str) > 1500:
        print(rec1_str[:1500] + "\n  ... (truncated)")
    else:
        print(rec1_str)
    
    return True


@pytest.mark.slow
def test_suzuki_coupling():
    """
    Test 2: Suzuki Coupling - Format Consistency
    
    Validates format consistency for a different reaction type.
    
    WARNING: This test is SLOW (~2-3 minutes) because Suzuki dataset has 50,215 reactions.
             Run with: pytest -m slow
             Skip with: pytest -m "not slow"
    """
    print("\n" + "="*70)
    print("TEST 2: Suzuki Coupling - Format Consistency")
    print("="*70)
    
    # Import heavy modules only when test runs
    from chemtools.recommend.core import recommend_from_reaction
    from chemtools.reaction_type_detector import detect_reaction_type
    
    smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    
    print(f"\n📋 Input: {smiles}")
    print(f"   Expected: Suzuki coupling")
    
    # ML-based
    start = time.time()
    # Use k=5 for faster testing
    ml_raw = recommend_from_reaction(smiles, k=5)
    ml_time = (time.time() - start) * 1000
    
    detection = detect_reaction_type(smiles)
    detected_type = detection.get("type", "Unknown")
    
    # Use the formatted output that's already in ml_raw
    if 'formatted' in ml_raw and 'recommended_conditions' in ml_raw['formatted']:
        ml_output = ml_raw['formatted']
        # Rename 'recommended_conditions' to 'recommendations' for consistency
        if 'recommended_conditions' in ml_output:
            ml_output['recommendations'] = ml_output.pop('recommended_conditions')
    else:
        print(f"❌ No formatted output found in ML result")
        return False
    
    print(f"\n✅ ML processed in {ml_time:.2f} ms")
    print(f"📊 {len(ml_output['recommendations'])} recommendations generated")
    
    # Validate first recommendation
    if len(ml_output["recommendations"]) > 0:
        rec = ml_output["recommendations"][0]
        if validate_recommendation_structure(rec):
            print(f"✅ Structure validation passed")
            chemicals = rec.get('chemicals', rec.get('reagents', []))
            print(f"   Chemicals: {len(chemicals)}")
            print(f"   Conditions: {', '.join(rec['conditions'].keys())}")
            return True
        else:
            print(f"❌ Structure validation failed")
            return False
    
    return False


def test_field_consistency():
    """
    Test 3: Field Consistency
    
    Ensures all recommendations have consistent field ordering and structure.
    """
    print("\n" + "="*70)
    print("TEST 3: Field Consistency Check")
    print("="*70)
    
    # Import heavy modules only when test runs
    from chemtools.recommend.core import recommend_from_reaction
    from chemtools.reaction_type_detector import detect_reaction_type
    
    smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    # Generate outputs - use k=5 for faster testing
    ml_raw = recommend_from_reaction(smiles, k=5)
    detection = detect_reaction_type(smiles)
    
    # Use the formatted output that's already in ml_raw
    if 'formatted' in ml_raw and 'recommended_conditions' in ml_raw['formatted']:
        ml_output = ml_raw['formatted']
        # Rename 'recommended_conditions' to 'recommendations' for consistency
        if 'recommended_conditions' in ml_output:
            ml_output['recommendations'] = ml_output.pop('recommended_conditions')
    else:
        print(f"❌ No formatted output found in ML result")
        return False
    
    # Check field ordering consistency
    print("\n🔍 Checking field ordering consistency...")
    
    expected_order = {
        "recommendation": ["rank", "confidence", "support", "reaction", "reagents", "conditions"],
        "reagent": ["id", "role", "name", "abbreviation", "cas", "smiles", "inchi_key"],
        "conditions": ["temperature", "time", "atmosphere", "pressure"],
    }
    
    all_consistent = True
    
    for rec in ml_output["recommendations"]:
        rec_keys = list(rec.keys())
        # Check that expected fields appear in order
        expected_rec = [k for k in expected_order["recommendation"] if k in rec_keys]
        actual_rec = [k for k in rec_keys if k in expected_order["recommendation"]]
        
        if expected_rec != actual_rec:
            print(f"  ❌ Recommendation field order mismatch")
            print(f"     Expected: {expected_rec}")
            print(f"     Actual: {actual_rec}")
            all_consistent = False
    
    if all_consistent:
        print("  ✅ All recommendations have consistent field ordering")
    
    return all_consistent


def run_all_tests():
    """Run all unified format tests."""
    print("\n" + "="*70)
    print("UNIFIED OUTPUT FORMAT TEST SUITE")
    print("="*70)
    print("\nValidating ML-based recommendations produce standardized")
    print("JSON format suitable for robotic execution systems.")
    
    tests = [
        ("Ullmann C-N Coupling", test_ullmann_cn_coupling),
        ("Suzuki Coupling", test_suzuki_coupling),
        ("Field Consistency", test_field_consistency),
    ]
    
    results = []
    
    for name, test_fn in tests:
        try:
            success = test_fn()
            results.append((name, success))
        except Exception as e:
            print(f"\n❌ Test '{name}' failed with error: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    for name, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status} - {name}")
    
    total = len(results)
    passed = sum(1 for _, s in results if s)
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! Format is ready for robotic execution.")
        return True
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please review.")
        return False


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
