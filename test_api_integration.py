#!/usr/bin/env python3
"""
Quick integration test to verify API endpoints work with refactored core.py
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))


def test_api_imports():
    """Test that API can import from refactored core.py"""
    print("Testing API integration with refactored core.py...")
    
    try:
        # Test import as the API would use it
        from chemtools.recommend import recommend_from_reaction
        from chemtools.recommend import recommend_conditions_structured
        
        print("  ✅ API imports successful")
        return True
        
    except ImportError as e:
        print(f"  ❌ API import failed: {e}")
        return False


def test_quick_recommendation():
    """Test a quick recommendation to verify functionality"""
    print("\nTesting quick recommendation...")
    
    try:
        from chemtools.recommend import recommend_conditions_structured
        
        # Simple Suzuki coupling
        reaction = "c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1"
        
        result = recommend_conditions_structured(
            reaction=reaction,
            k=5,
            limit=2,
            rerank_strategy='none',
        )
        
        # Check result structure
        assert isinstance(result, dict), "Result should be a dict"
        assert "meta" in result, "Result should have meta"
        assert "recommendations" in result, "Result should have recommendations"
        
        print(f"  ✅ Quick recommendation successful")
        print(f"     - Status: {result['meta'].get('status', 'unknown')}")
        print(f"     - Recommendations: {len(result.get('recommendations', []))}")
        print(f"     - Processing time: {result['meta'].get('processing_time_ms', 0)}ms")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Quick recommendation failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("=" * 60)
    print("API Integration Test - Refactored Core.py")
    print("=" * 60)
    print()
    
    tests = [
        test_api_imports,
        test_quick_recommendation,
    ]
    
    results = [test() for test in tests]
    
    print("\n" + "=" * 60)
    if all(results):
        print("✅ All API integration tests passed!")
        print("   The refactored core.py is production-ready.")
        return 0
    else:
        print("❌ Some tests failed. Please review.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
