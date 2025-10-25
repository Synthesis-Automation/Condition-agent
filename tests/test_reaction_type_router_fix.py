"""
Test script to verify the reaction type router fix.

Tests that C-N coupling reactions are correctly routed to cn_coupling_Cu_db.json
instead of incorrectly using Suzuki_db.json.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.local_recommendation_cli import local_rule_based_match
from scripts.recommendation_cli_utils import REACTION_TYPE_CHOICES
import pytest
from chemtools.router import detect_family_from_reaction

SCDB_DIR = Path("data/rule_db")

def _ensure_rule_db(filenames):
    missing = [SCDB_DIR / name for name in filenames if not (SCDB_DIR / name).exists()]
    if missing:
        pytest.skip(f"Rule DB not available: {missing}")



def test_reaction_type_choices():
    """Verify that deprecated catalyst-specific options are removed."""
    print("=" * 70)
    print("TEST 1: Reaction Type Choices")
    print("=" * 70)
    
    # Should have unified cn_coupling, not separate Cu/Pd/Ni
    choices = [choice[1] for choice in REACTION_TYPE_CHOICES if choice[1]]
    
    print(f"Available choices: {choices}")
    
    # Check for unified option
    assert "cn_coupling" in choices, "cn_coupling should be available"
    
    # Check deprecated options are removed
    deprecated = ["C_N_Coupling_Cu", "C_N_Coupling_Pd", "C_N_Coupling_Ni"]
    for dep in deprecated:
        assert dep not in choices, f"Deprecated option {dep} should be removed"
    
    print("✅ PASSED: Reaction type choices are correct")
    print()


def test_router_detection():
    """Verify that router correctly detects C-N coupling."""
    print("=" * 70)
    print("TEST 2: Router Detection")
    print("=" * 70)
    
    # User's example: aryl bromide + benzylamine
    rsmi = "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1"
    
    result = detect_family_from_reaction(rsmi, use_rxn_insight=False)
    
    family = result.get("family")
    confidence = result.get("confidence")
    
    print(f"Reaction: {rsmi}")
    print(f"Detected Family: {family}")
    print(f"Confidence: {confidence}")
    
    assert family == "cn_coupling", f"Expected cn_coupling, got {family}"
    assert confidence >= 0.85, f"Expected confidence >= 0.85, got {confidence}"
    
    print("✅ PASSED: Router correctly detects C-N coupling")
    print()


def test_database_selection():
    """Verify that correct database is selected for C-N coupling."""
    print("=" * 70)
    print("TEST 3: Database Selection")
    print("=" * 70)
    
    rsmi = "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1"
    reaction_type = "cn_coupling"
    
    _ensure_rule_db(["C_N_Coupling_Cu_db.json"])
    result = local_rule_based_match(rsmi, None, reaction_type)
    
    # Check no errors
    assert "error" not in result, f"Got error: {result.get('error')}"
    
    detection = result.get("detection", {})
    family = detection.get("family")
    
    print(f"Reaction Type: {reaction_type}")
    print(f"Detected Family: {family}")
    
    # The family might be ullmann_cn or cn_coupling (both are valid for C-N coupling)
    valid_families = ["cn_coupling", "ullmann_cn", "buchwald_hartwig_c_n"]
    assert family in valid_families, f"Expected C-N coupling variant, got {family}"
    
    # Check recommendations exist
    recommendations = result.get("recommended_conditions", [])
    assert len(recommendations) > 0, "Should have at least one recommendation"
    
    top = recommendations[0]
    source = top.get("source", {})
    entry_name = source.get("entry_name", "")
    
    print(f"Top Match: {entry_name}")
    print(f"Recommendations: {len(recommendations)} found")
    
    # The entry should be related to C-N coupling (ArBr + amine)
    # Not Suzuki-related (ArBr + boronic acid)
    assert "amine" in entry_name.lower() or "ArBr" in entry_name, \
        f"Expected C-N coupling entry, got: {entry_name}"
    
    print("✅ PASSED: Correct database selected and recommendations found")
    print()


def test_suzuki_still_works():
    """Verify that Suzuki detection still works correctly."""
    print("=" * 70)
    print("TEST 4: Suzuki Detection (Regression Test)")
    print("=" * 70)
    
    # Suzuki reaction: aryl bromide + boronic acid
    rsmi = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    
    result_auto = detect_family_from_reaction(rsmi, use_rxn_insight=False)
    family_auto = result_auto.get("family")
    
    print(f"Reaction: {rsmi}")
    print(f"Auto-detected Family: {family_auto}")
    
    # Should detect as Suzuki
    assert family_auto == "suzuki_miyaura", f"Expected Suzuki, got {family_auto}"
    
    # Test rule-based matching
    _ensure_rule_db(["Suzuki_db.json"])
    result_rule = local_rule_based_match(rsmi, None, "suzuki_miyaura")
    
    assert "error" not in result_rule, f"Got error: {result_rule.get('error')}"
    
    recommendations = result_rule.get("recommended_conditions", [])
    assert len(recommendations) > 0, "Should have Suzuki recommendations"
    
    print(f"Recommendations: {len(recommendations)} found")
    print("✅ PASSED: Suzuki detection still works")
    print()


def test_comparison_before_after():
    """Compare behavior with different reaction types."""
    print("=" * 70)
    print("TEST 5: Comparison - C-N Coupling vs Suzuki")
    print("=" * 70)
    
    cn_reaction = "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1"
    suzuki_reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    
    # Test C-N coupling
    _ensure_rule_db(["C_N_Coupling_Cu_db.json"])
    cn_result = local_rule_based_match(cn_reaction, None, "cn_coupling")
    cn_family = cn_result.get("detection", {}).get("family", "Unknown")
    cn_matches = len(cn_result.get("recommended_conditions", []))
    
    # Test Suzuki
    _ensure_rule_db(["Suzuki_db.json"])
    suzuki_result = local_rule_based_match(suzuki_reaction, None, "suzuki_miyaura")
    suzuki_family = suzuki_result.get("detection", {}).get("family", "Unknown")
    suzuki_matches = len(suzuki_result.get("recommended_conditions", []))
    
    print(f"C-N Coupling:")
    print(f"  Reaction: ArBr + Amine")
    print(f"  Detected: {cn_family}")
    print(f"  Matches: {cn_matches}")
    
    print(f"\nSuzuki:")
    print(f"  Reaction: ArBr + Boronic Acid")
    print(f"  Detected: {suzuki_family}")
    print(f"  Matches: {suzuki_matches}")
    
    # They should be different families
    assert cn_family != suzuki_family, \
        "C-N coupling should not be detected as Suzuki"
    
    print("\n✅ PASSED: Different reactions routed to different databases")
    print()


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("REACTION TYPE ROUTER FIX - VERIFICATION TESTS")
    print("=" * 70)
    print()
    
    try:
        test_reaction_type_choices()
        test_router_detection()
        test_database_selection()
        test_suzuki_still_works()
        test_comparison_before_after()
        
        print("=" * 70)
        print("ALL TESTS PASSED ✅")
        print("=" * 70)
        print("\nThe reaction type router fix is working correctly:")
        print("  ✓ Deprecated catalyst-specific options removed")
        print("  ✓ Router detects C-N coupling correctly")
        print("  ✓ Correct database selected for C-N coupling")
        print("  ✓ Suzuki detection still works")
        print("  ✓ Different reactions routed to different databases")
        print()
        
    except AssertionError as e:
        print("\n" + "=" * 70)
        print("TEST FAILED ❌")
        print("=" * 70)
        print(f"\nError: {e}")
        print()
        sys.exit(1)
    except Exception as e:
        print("\n" + "=" * 70)
        print("UNEXPECTED ERROR ❌")
        print("=" * 70)
        print(f"\nError: {type(e).__name__}: {e}")
        print()
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
