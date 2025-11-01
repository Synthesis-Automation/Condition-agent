"""
Test ChemTools wrapper functions directly (without LLM).

This script tests the tool wrappers by calling them directly,
bypassing the LangChain agent. Useful for verification without
needing API keys.

Usage:
    python -m lang_chain.test_tools
"""

import json
import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chem_assistant.chemtools_wrapper import (
    normalize_smiles_tool,
    normalize_reaction_tool,
    detect_reaction_family_tool,
    classify_reactant_tool,
    get_functional_groups_tool,
    recommend_conditions_tool,
    enhanced_cross_family_recommend_tool,
    search_precedents_tool,
    find_reagent_tool,
    list_supported_cores_tool,
    add_reagent_tool,
)


def _format_output(payload: object, limit: int = 400) -> str:
    """Pretty-print helpers for consistent console output."""
    try:
        text = json.dumps(payload, indent=2, default=str)
    except TypeError:
        text = str(payload)
    if len(text) > limit:
        return text[:limit] + "..."
    return text


def _assert_success(result: dict, message: str) -> None:
    """Validate standard tool success response."""
    assert isinstance(result, dict), f"{message}: tool did not return a dict"
    assert result.get("success"), f"{message}: {result}"


def test_normalize_smiles():
    """Test SMILES normalization."""
    print("\n" + "="*70)
    print("Test 1: Normalize SMILES")
    print("="*70)
    
    result = normalize_smiles_tool.invoke({"smiles": "c1ccccc1"})
    print(f"Input: c1ccccc1")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "SMILES normalization failed")
    assert result.get("smiles_norm") == "c1ccccc1", "Missing normalized SMILES"
    print("PASSED")


def test_normalize_reaction():
    """Test reaction normalization."""
    print("\n" + "="*70)
    print("Test 2: Normalize Reaction")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = normalize_reaction_tool.invoke({"reaction_smiles": reaction})
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Reaction normalization failed")
    assert result.get("normalized"), "Normalized reaction missing"
    print("PASSED")


def test_detect_reaction_family():
    """Test reaction family detection."""
    print("\n" + "="*70)
    print("Test 3: Detect Reaction Family")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = detect_reaction_family_tool.invoke({"reaction_smiles": reaction})
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Family detection failed")
    assert result.get("family"), "Family label missing"
    print("PASSED")


def test_classify_reactant():
    """Test reactant classification."""
    print("\n" + "="*70)
    print("Test 4: Classify Reactant")
    print("="*70)
    
    result = classify_reactant_tool.invoke({"smiles": "Brc1ccccc1"})
    print(f"Input: Brc1ccccc1")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Reactant classification failed")
    assert "classification" in result or "classes" in result, "Classification metadata missing"
    print("PASSED")


def test_functional_groups():
    """Test functional group detection."""
    print("\n" + "="*70)
    print("Test 5: Detect Functional Groups")
    print("="*70)
    
    result = get_functional_groups_tool.invoke({"smiles": "CCO"})
    print(f"Input: CCO")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Functional group detection failed")
    assert "alcohol" in result, "Expected alcohol functional group metadata"
    print("PASSED")


def test_find_reagent():
    """Test reagent lookup."""
    print("\n" + "="*70)
    print("Test 6: Find Reagent")
    print("="*70)
    
    result = find_reagent_tool.invoke({"query": "Cs2CO3"})
    print(f"Input: Cs2CO3")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Reagent lookup failed")
    assert result.get("cas"), "CAS metadata missing"
    print("PASSED")


def test_recommend_conditions():
    """Test condition recommendations."""
    print("\n" + "="*70)
    print("Test 7: Recommend Conditions")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = recommend_conditions_tool.invoke({
        "reaction_smiles": reaction,
        "k": 10,
        "max_variants": 2
    })
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Recommendation failed")
    assert "recommendation" in result, "Recommendation payload missing"
    assert "cache_hit" in result and "timing_ms" in result
    print("PASSED")


def test_enhanced_cross_family():
    """Test enhanced cross-family recommendation workflow."""
    print("\n" + "="*70)
    print("Test 8: Enhanced Cross-Family Recommendations")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = enhanced_cross_family_recommend_tool.invoke({
        "reaction_smiles": reaction,
        "k": 10,
    })
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Enhanced cross-family recommendation failed")
    metrics = result.get("cross_family_metrics", {})
    assert metrics, "Cross-family metrics missing"
    assert result.get("insights"), "Insights should be returned for LLM guidance"
    assert metrics.get("total_recommendations", 0) >= 1, "No cross-family precedents reported"
    print("PASSED")


def test_search_precedents():
    """Test precedent search."""
    print("\n" + "="*70)
    print("Test 9: Search Precedents")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = search_precedents_tool.invoke({
        "reaction_smiles": reaction,
        "k": 5
    })
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Precedent search failed")
    assert result.get("precedents"), "No precedents returned"
    print("PASSED")


def test_recommend_with_constraints():
    """Ensure constraint-aware parameters are accepted."""
    print("\n" + "="*70)
    print("Test 10: Recommend With Constraints")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = recommend_conditions_tool.invoke({
        "reaction_smiles": reaction,
        "constraint_text": "Pd-free, prefer copper",
        "allow_metals": ["Cu", "Ni"],
        "constraint_rules": {"no_chlorinated": True},
    })
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Constraint-aware recommendation failed")
    assert "constraint_summary" in result or "constraint_notes" in result, "Constraint metadata missing"
    assert "cache_hit" in result
    print("PASSED")


def test_list_supported_cores():
    """Test the core listing helper tool."""
    print("\n" + "="*70)
    print("Test 11: List Supported Cores")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = list_supported_cores_tool.invoke({
        "reaction_smiles": reaction,
        "k": 5
    })
    print(f"Input: {reaction}")
    print(f"Output: {_format_output(result)}")
    _assert_success(result, "Core listing failed")
    assert "core_candidates" in result, "Core candidates missing"
    assert "cache_hit" in result and "timing_ms" in result
    print("PASSED")

def test_add_reagent_tool_dry_run():
    """Ensure reagent addition tool responds without crashing."""
    print("\n" + "="*70)
    print("Test 12: Add Reagent (Dry Run)")
    print("="*70)
    
    payload = {
        "cas": "50-00-0",
        "name": "Formaldehyde",
        "role": "other_reagent",
        "allow_default_family": True,
        "dry_run": True,
        "auto_resolve": False,
    }
    result = add_reagent_tool.invoke(payload)
    print(f"Payload: {payload}")
    print(f"Output: {_format_output(result)}")
    if result.get("success"):
        assert result.get("status") in {"dry_run", "exists", "written", None}, "Unexpected reagent status"
    else:
        print("Warning: add_reagent_tool reported an error (acceptable in dry-run environments).")
    print("PASSED")


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("ChemTools LangChain Wrapper - Direct Tool Tests")
    print("="*70)
    print("\nTesting tools without LLM agent (no API key needed)")
    
    tests = [
        test_normalize_smiles,
        test_normalize_reaction,
        test_detect_reaction_family,
        test_classify_reactant,
        test_functional_groups,
        test_find_reagent,
        test_recommend_conditions,
        test_enhanced_cross_family,
        test_search_precedents,
        test_recommend_with_constraints,
        test_list_supported_cores,
        test_add_reagent_tool_dry_run,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"FAILED: {e}")
            failed += 1
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*70)
    
    if failed == 0:
        print("\nAll tests passed! The wrapper is working correctly.")
        print("\nNext steps:")
        print("  1. Set API key: $env:OPENAI_API_KEY = 'sk-your-key-here'")
        print("  2. Run CLI: python -m lang_chain.chemtools_cli")
        print("  3. Or run examples: python lang_chain/example_usage.py")
    else:
        print(f"\n {failed} test(s) failed. Check the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
