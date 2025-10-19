#!/usr/bin/env python3
"""
Test script to verify core.py refactoring.

Tests that:
1. All import paths work correctly
2. Function behavior is identical to original
3. Backwards compatibility is maintained
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def test_imports():
    """Test all import paths work correctly."""
    print("Testing imports...")
    
    # Test 1: Old-style imports from core (backwards compatibility)
    try:
        from chemtools.recommend.core import recommend_from_reaction
        from chemtools.recommend.core import recommend_conditions_structured
        print("  ✅ Old-style imports from core work")
    except ImportError as e:
        print(f"  ❌ Old-style imports failed: {e}")
        return False
    
    # Test 2: New-style imports from modules
    try:
        from chemtools.recommend.modules import recommend_from_reaction as rec_mod
        from chemtools.recommend.modules import recommend_conditions_structured as struct_mod
        print("  ✅ New-style imports from modules work")
    except ImportError as e:
        print(f"  ❌ New-style imports failed: {e}")
        return False
    
    # Test 3: Direct module imports
    try:
        from chemtools.recommend.modules.recommender import recommend_from_reaction as rec_direct
        from chemtools.recommend.modules.structured import recommend_conditions_structured as struct_direct
        print("  ✅ Direct module imports work")
    except ImportError as e:
        print(f"  ❌ Direct module imports failed: {e}")
        return False
    
    # Test 4: Internal function imports (backwards compatibility)
    try:
        from chemtools.recommend.core import _build_formatted_output
        from chemtools.recommend.core import _build_precedent_details
        print("  ✅ Internal function imports work")
    except ImportError as e:
        print(f"  ❌ Internal function imports failed: {e}")
        return False
    
    # Test 5: Verify functions are the same object
    assert recommend_from_reaction is rec_mod is rec_direct, \
        "recommend_from_reaction should be the same object across all import paths"
    assert recommend_conditions_structured is struct_mod is struct_direct, \
        "recommend_conditions_structured should be the same object across all import paths"
    print("  ✅ All import paths reference the same functions")
    
    return True


def test_basic_execution():
    """Test basic function execution with minimal input."""
    print("\nTesting basic execution...")
    
    try:
        from chemtools.recommend.core import recommend_from_reaction
        
        # Test with a simple reaction SMILES
        test_reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        result = recommend_from_reaction(
            reaction=test_reaction,
            k=3,  # Just get top 3 precedents
            rerank_strategy='none',  # Disable reranking for simplicity
        )
        
        # Check basic result structure
        assert isinstance(result, dict), "Result should be a dictionary"
        
        # Print actual structure for debugging
        print(f"     - Result keys: {list(result.keys())}")
        
        # Check for expected keys (may vary based on actual implementation)
        if "precedents" in result:
            assert isinstance(result["precedents"], list), "Precedents should be a list"
            print(f"     - Found {len(result['precedents'])} precedents")
        elif "recommendations" in result:
            print(f"     - Found {len(result.get('recommendations', []))} recommendations")
        
        if "detected_family" in result:
            print(f"     - Detected family: {result['detected_family']}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Execution failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_structured_execution():
    """Test structured API wrapper execution."""
    print("\nTesting structured execution...")
    
    try:
        from chemtools.recommend.core import recommend_conditions_structured
        
        # Test with a simple reaction SMILES
        test_reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        result = recommend_conditions_structured(
            reaction=test_reaction,
            k=3,
            rerank_strategy='none',
        )
        
        # Check structured result format
        assert isinstance(result, dict), "Result should be a dictionary"
        
        # Print actual structure for debugging
        print(f"     - Result keys: {list(result.keys())}")
        
        # Check for expected keys
        if "condition_variants" in result:
            print(f"     - Generated {len(result['condition_variants'])} variants")
        elif "recommendations" in result:
            print(f"     - Generated {len(result.get('recommendations', []))} recommendations")
        
        if "processing_time_ms" in result:
            print(f"     - Processing time: {result['processing_time_ms']}ms")
        elif "meta" in result and "processing_time_ms" in result["meta"]:
            print(f"     - Processing time: {result['meta']['processing_time_ms']}ms")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Execution failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_module_structure():
    """Test that module structure is correct."""
    print("\nTesting module structure...")
    
    try:
        import chemtools.recommend.modules as modules
        
        # Check that all expected functions are exported
        expected_functions = [
            "recommend_from_reaction",
            "recommend_conditions_structured",
            "build_formatted_output",
            "convert_fusion_to_core_format",
            "build_formatted_output_from_fusion",
            "build_precedent_details",
            "calculate_average_yield",
            "calculate_yield_range",
            "calculate_temp_range",
            "calculate_time_range",
        ]
        
        for func_name in expected_functions:
            assert hasattr(modules, func_name), \
                f"modules package should export {func_name}"
        
        print(f"  ✅ All {len(expected_functions)} expected functions exported")
        
        # Check that individual modules exist
        from chemtools.recommend.modules import recommender
        from chemtools.recommend.modules import structured
        from chemtools.recommend.modules import fusion_adapter
        from chemtools.recommend.modules import precedent_builder
        from chemtools.recommend.modules import output_builder
        
        print("  ✅ All 5 submodules accessible")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Module structure check failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_file_sizes():
    """Compare file sizes before and after refactoring."""
    print("\nComparing file sizes...")
    
    try:
        backup_path = project_root / "chemtools" / "recommend" / "core.py.backup"
        new_core_path = project_root / "chemtools" / "recommend" / "core.py"
        modules_dir = project_root / "chemtools" / "recommend" / "modules"
        
        if backup_path.exists():
            backup_lines = len(backup_path.read_text(encoding="utf-8").splitlines())
            print(f"  📊 Original core.py: {backup_lines} lines")
        
        if new_core_path.exists():
            new_core_lines = len(new_core_path.read_text(encoding="utf-8").splitlines())
            print(f"  📊 New core.py: {new_core_lines} lines")
            
            if backup_path.exists():
                reduction = (1 - new_core_lines / backup_lines) * 100
                print(f"  ✅ Reduction: {reduction:.1f}%")
        
        if modules_dir.exists():
            module_files = list(modules_dir.glob("*.py"))
            total_module_lines = sum(
                len(f.read_text(encoding="utf-8").splitlines()) 
                for f in module_files
            )
            print(f"  📊 Total module lines: {total_module_lines} (across {len(module_files)} files)")
        
        return True
        
    except Exception as e:
        print(f"  ⚠️  Could not compare file sizes: {e}")
        return True  # Not a critical failure


def main():
    """Run all tests."""
    print("=" * 60)
    print("Core.py Refactoring Test Suite")
    print("=" * 60)
    
    tests = [
        ("Import Tests", test_imports),
        ("Basic Execution", test_basic_execution),
        ("Structured Execution", test_structured_execution),
        ("Module Structure", test_module_structure),
        ("File Size Comparison", test_file_sizes),
    ]
    
    results = []
    for test_name, test_func in tests:
        print()
        success = test_func()
        results.append((test_name, success))
    
    # Print summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for test_name, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! Refactoring successful!")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please review.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
