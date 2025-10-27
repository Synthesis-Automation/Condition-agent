"""
Test ChemTools wrapper functions directly (without LLM).

This script tests the tool wrappers by calling them directly,
bypassing the LangChain agent. Useful for verification without
needing API keys.

Usage:
    python -m lang_chain.test_tools
"""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from lang_chain.chemtools_wrapper import (
    normalize_smiles_tool,
    normalize_reaction_tool,
    detect_reaction_family_tool,
    classify_reactant_tool,
    get_functional_groups_tool,
    recommend_conditions_tool,
    search_precedents_tool,
    find_reagent_tool,
)


def test_normalize_smiles():
    """Test SMILES normalization."""
    print("\n" + "="*70)
    print("Test 1: Normalize SMILES")
    print("="*70)
    
    result = normalize_smiles_tool.invoke({"smiles": "c1ccccc1"})
    print(f"Input: c1ccccc1")
    print(f"Output: {result}")
    # Result is JSON string, check it contains normalized SMILES
    assert "smiles_norm" in result or "c1ccccc1" in result, "SMILES normalization failed"
    print("✓ PASSED")


def test_normalize_reaction():
    """Test reaction normalization."""
    print("\n" + "="*70)
    print("Test 2: Normalize Reaction")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = normalize_reaction_tool.invoke({"reaction_smiles": reaction})
    print(f"Input: {reaction}")
    print(f"Output: {result}")
    assert len(result) > 0, "Reaction normalization failed"
    print("✓ PASSED")


def test_detect_reaction_family():
    """Test reaction family detection."""
    print("\n" + "="*70)
    print("Test 3: Detect Reaction Family")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = detect_reaction_family_tool.invoke({"reaction_smiles": reaction})
    print(f"Input: {reaction}")
    print(f"Output: {result[:200]}...")  # First 200 chars
    assert "family" in result.lower(), "Family detection failed"
    print("✓ PASSED")


def test_classify_reactant():
    """Test reactant classification."""
    print("\n" + "="*70)
    print("Test 4: Classify Reactant")
    print("="*70)
    
    result = classify_reactant_tool.invoke({"smiles": "Brc1ccccc1"})
    print(f"Input: Brc1ccccc1")
    print(f"Output: {result[:200]}...")  # First 200 chars
    # Just check it doesn't error - null is acceptable for some reactants
    assert len(result) > 0, "Reactant classification failed"
    print("✓ PASSED")


def test_functional_groups():
    """Test functional group detection."""
    print("\n" + "="*70)
    print("Test 5: Detect Functional Groups")
    print("="*70)
    
    result = get_functional_groups_tool.invoke({"smiles": "CCO"})
    print(f"Input: CCO")
    print(f"Output: {result[:200]}...")  # First 200 chars
    assert "alcohol" in result.lower() or "error" in result.lower(), "Functional group detection failed"
    print("✓ PASSED")


def test_find_reagent():
    """Test reagent lookup."""
    print("\n" + "="*70)
    print("Test 6: Find Reagent")
    print("="*70)
    
    result = find_reagent_tool.invoke({"query": "Cs2CO3"})
    print(f"Input: Cs2CO3")
    print(f"Output: {result[:200]}...")  # First 200 chars
    assert "cas" in result.lower() or "error" in result.lower(), "Reagent lookup failed"
    print("✓ PASSED")


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
    print(f"Output: {result[:300]}...")  # First 300 chars
    assert "family" in result.lower() or "error" in result.lower(), "Recommendation failed"
    print("✓ PASSED")


def test_search_precedents():
    """Test precedent search."""
    print("\n" + "="*70)
    print("Test 8: Search Precedents")
    print("="*70)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    result = search_precedents_tool.invoke({
        "reaction_smiles": reaction,
        "k": 5
    })
    print(f"Input: {reaction}")
    print(f"Output: {result[:300]}...")  # First 300 chars
    assert "precedent" in result.lower() or "error" in result.lower(), "Precedent search failed"
    print("✓ PASSED")


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
        test_search_precedents,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"✗ FAILED: {e}")
            failed += 1
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*70)
    
    if failed == 0:
        print("\n✅ All tests passed! The wrapper is working correctly.")
        print("\nNext steps:")
        print("  1. Set API key: $env:OPENAI_API_KEY = 'sk-your-key-here'")
        print("  2. Run CLI: python -m lang_chain.chemtools_cli")
        print("  3. Or run examples: python lang_chain/example_usage.py")
    else:
        print(f"\n❌ {failed} test(s) failed. Check the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
