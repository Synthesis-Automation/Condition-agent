#!/usr/bin/env python
"""Run only Test 8 for catalytic system generation analysis."""
import sys
from pathlib import Path

parent_dir = Path(__file__).parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

# Add tests directory to path
tests_dir = parent_dir / "tests"
if str(tests_dir) not in sys.path:
    sys.path.insert(0, str(tests_dir))

# Import the test
import test_analytics_module

if __name__ == "__main__":
    try:
        success = test_analytics_module.test_8_catalytic_systems()
        print(f"\n{'='*80}")
        if success:
            print("✅ Test 8 completed successfully!")
        else:
            print("❌ Test 8 failed!")
        print('='*80)
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ Error running test: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
