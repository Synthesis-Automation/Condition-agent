"""
Fast version of unified output format test - skips slow Suzuki test.

This version only tests:
- Test 1: Ullmann C-N Coupling (~1 second)
- Test 3: Field Consistency (~1 second)

Total run time: ~2-3 seconds instead of ~3 minutes

For full test suite including Suzuki, run:
    python tests/test_unified_output_format.py
"""

import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Add tests directory to path
tests_dir = Path(__file__).parent
if str(tests_dir) not in sys.path:
    sys.path.insert(0, str(tests_dir))

# Import the test functions
from test_unified_output_format import (
    test_ullmann_cn_coupling,
    test_field_consistency,
    run_all_tests
)


def run_fast_tests():
    """Run only the fast tests, skip Suzuki."""
    print("=" * 70)
    print("FAST UNIFIED OUTPUT FORMAT TEST SUITE")
    print("=" * 70)
    print()
    print("Running fast tests only (skipping slow Suzuki test)")
    print()
    
    tests = [
        ("Ullmann C-N Coupling", test_ullmann_cn_coupling),
        ("Field Consistency", test_field_consistency),
    ]
    
    results = {}
    for name, test_fn in tests:
        try:
            success = test_fn()
            results[name] = "PASS" if success else "FAIL"
        except Exception as e:
            print(f"\n❌ Test '{name}' failed with error: {e}")
            import traceback
            traceback.print_exc()
            results[name] = "FAIL"
    
    # Summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    for name, status in results.items():
        icon = "✅" if status == "PASS" else "❌"
        print(f"{icon} {status} - {name}")
    
    passed = sum(1 for s in results.values() if s == "PASS")
    total = len(results)
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All fast tests passed! Format is ready for robotic execution.")
        print("\n💡 To run full test suite including Suzuki (~3 min):")
        print("   python tests/test_unified_output_format.py")
        return True
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please review.")
        return False


if __name__ == "__main__":
    success = run_fast_tests()
    sys.exit(0 if success else 1)
