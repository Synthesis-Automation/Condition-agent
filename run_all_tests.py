#!/usr/bin/env python3
"""
Simple test runner for all core.py refactoring tests.
Runs all test suites and provides a summary.
"""

import subprocess
import sys
from pathlib import Path

def run_command(cmd, description):
    """Run a command and return success status"""
    print(f"\n{'=' * 70}")
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    print('=' * 70)
    
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        # Print output
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        success = result.returncode == 0
        status = "[PASS]" if success else "[FAIL]"
        print(f"\n{status} Exit code: {result.returncode}")
        
        return success
        
    except subprocess.TimeoutExpired:
        print("[FAIL] Command timed out after 120 seconds")
        return False
    except Exception as e:
        print(f"[FAIL] Error running command: {e}")
        return False


def main():
    """Run all test suites"""
    print("=" * 70)
    print("CORE.PY REFACTORING - COMPLETE TEST SUITE")
    print("=" * 70)
    
    tests = [
        ("python test_core_refactoring.py", "Core Refactoring Tests"),
        ("pytest tests/test_core_modules.py -v --tb=short", "Unit Tests (pytest)"),
        ("python test_api_integration.py", "API Integration Tests"),
    ]
    
    results = []
    for cmd, description in tests:
        success = run_command(cmd, description)
        results.append((description, success))
    
    # Print summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for description, success in results:
        status = "[PASS]" if success else "[FAIL]"
        print(f"{status} {description}")
    
    print(f"\nTotal: {passed}/{total} test suites passed")
    
    if passed == total:
        print("\n[SUCCESS] All test suites passed!")
        print("The refactored core.py is production-ready.")
        return 0
    else:
        print(f"\n[WARNING] {total - passed} test suite(s) failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
